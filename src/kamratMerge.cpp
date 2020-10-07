#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>
#include <stdexcept>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "utils/seq_coding.hpp"
#include "data_struct/count_tab.hpp"
#include "data_struct/contig_elem.hpp"
#include "data_struct/merge_knot.hpp"
#include "run_info_parser/merge.hpp"
#include "run_info_parser/utils.hpp"

const size_t EstimateAvgCount(std::vector<float> &count_avg_vect,
                              const CountTab &kmer_count_tab,
                              const code2serial_t code2serial,
                              const std::string &seq,
                              const size_t k_len, const bool stranded, const size_t max_dist,
                              const std::string &interv_method, const float interv_thres,
                              std::ifstream &idx_file)
{
    size_t nb_kmer_found(1), seq_len = seq.size();
    if (seq_len < k_len) //shouldn't happen
    {
        throw std::domain_error("sequence size less than k-mer length");
    }
    if (seq_len == k_len)
    {
        auto iter = code2serial.find(Seq2Int(seq, k_len, stranded));
        if (iter == code2serial.cend())
        {
            throw std::domain_error("single-mer contig not found in k-mer count table");
        }
        kmer_count_tab.GetCountVect(count_avg_vect, iter->second, idx_file);
    }
    else
    {
        size_t pos_left(0), pos_right = pos_left;
        uint64_t uniqcode1 = Seq2Int(seq.substr(pos_left, k_len), k_len, stranded),
                 uniqcode2 = uniqcode1;
        auto iter1 = code2serial.find(uniqcode1);
        if (iter1 == code2serial.cend()) // shouldn't happen
        {
            throw std::domain_error("head k-mer not found in k-mer count table");
        }
        kmer_count_tab.GetCountVect(count_avg_vect, iter1->second, idx_file); // initialized with the head k-mer
        std::vector<float> count_vect_left(count_avg_vect);

        float count_dist(1);
        while (pos_left < seq_len - k_len) // pos_right <= seq_len - k_len
        {
            if (pos_right - pos_left > max_dist) // should not happen
            {
                throw std::domain_error("no adjacent k-mers found in sequence for given max distance");
            }
            ++pos_right;
            uniqcode2 = NextSeq(uniqcode2, k_len, seq[pos_right + k_len - 1], stranded);
            auto iter2 = code2serial.find(uniqcode2);
            if (iter2 == code2serial.cend())
            {
                count_dist = 1; // for debug
                continue;
            }
            count_dist = kmer_count_tab.AddCountVectIfCoherent(count_avg_vect, iter2->second, count_vect_left, interv_method, interv_thres, idx_file);
            if (count_dist < interv_thres) // merging is rejected when larger or equal
            {
                pos_left = pos_right;
                kmer_count_tab.GetCountVect(count_vect_left, iter2->second, idx_file);
                ++nb_kmer_found;
            }
        }
        if (count_dist >= interv_thres) // for debug, should not happen
        {
            throw std::domain_error("rear k-mer not found in k-mer count table");
        }
        for (unsigned int i(0); i < count_avg_vect.size(); ++i)
        {
            count_avg_vect[i] /= nb_kmer_found;
        }
    }
    return nb_kmer_found;
}

void ScanCountTable(CountTab &kmer_count_tab,
                    contigvect_t &contig_vect, code2serial_t &code2serial,
                    const size_t klen, const bool stranded,
                    const std::string &kmer_count_path,
                    const std::string &rep_colname,
                    const std::string &smp_info_path)
{
    std::ifstream kmer_count_file(kmer_count_path);
    if (!kmer_count_file.is_open())
    {
        throw std::domain_error("k-mer count file " + kmer_count_path + " was not found");
    }
    size_t pos = kmer_count_path.find_last_of(".");
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    if (pos != std::string::npos && kmer_count_path.substr(pos + 1) == "gz")
    {
        inbuf.push(boost::iostreams::gzip_decompressor());
    }
    inbuf.push(kmer_count_file);
    std::istream kmer_count_instream(&inbuf);

    //----- Dealing with Header Line for Constructing ColumnInfo Object -----//
    std::string line;
    std::getline(kmer_count_instream, line);
    kmer_count_tab.MakeSmpCond(smp_info_path);
    kmer_count_tab.MakeColumnInfo(line, rep_colname);
    //----- Dealing with Following k-mer Count Lines -----//
    std::string idx_path = kmer_count_tab.GetIndexPath();
    std::ofstream count_idx_file;
    if (!idx_path.empty())
    {
        count_idx_file.open(idx_path);
        if (!count_idx_file.is_open()) // to ensure the file is opened
        {
            throw std::domain_error("error open file: " + idx_path);
        }
    }
    for (size_t kmer_serial(0); std::getline(kmer_count_instream, line); ++kmer_serial)
    {
        float rep_val;
        bool has_rep_val;

        if (kmer_serial != kmer_count_tab.GetTableSize()) // for debug, should not happen
        {
            throw std::domain_error("kmer_count_tab not coherent with k-mer serial");
        }
        // parsing count vector => k-mer count table
        has_rep_val = kmer_count_tab.AddRowAsFields(rep_val, line, count_idx_file);
        if (!has_rep_val) // if no representative k-mer column
        {
            rep_val = kmer_serial;
        }
        // dealing with sequence => contig vector
        std::istringstream conv(line);
        std::string seq;
        conv >> seq; // first column as feature (string)
        if (seq.size() != klen)
        {
            throw std::domain_error("the given k-len parameter not coherent with input k-mer: " + seq);
        }
        uint64_t kmer_uniqcode = Seq2Int(seq, seq.size(), stranded);
        if (!code2serial.insert({kmer_uniqcode, kmer_serial}).second)
        {
            throw std::domain_error("duplicate input (newly inserted k-mer already exists in k-mer hash list)");
        }
        if (kmer_serial != contig_vect.size()) // for debug, should not happen
        {
            throw std::domain_error("contig vector not coherent with k-mer serial");
        }
        contig_vect.push_back(ContigElem(seq, rep_val, kmer_uniqcode, kmer_serial));
    }
    if (count_idx_file.is_open())
    {
        count_idx_file.close();
    }
    kmer_count_file.close();
}

void MakeOverlapKnotDict(fix2knot_t &hashed_mergeknot_list,
                         const contigvect_t &contig_vect,
                         const bool stranded,
                         const size_t n_overlap)
{
    for (size_t contig_serial(0); contig_serial < contig_vect.size(); ++contig_serial)
    {
        std::string seq = contig_vect[contig_serial].GetSeq();
        if (seq.size() <= n_overlap) // for debug, should not happen
        {
            throw std::domain_error("k-mer size less than or equal to k-length");
        }
        uint64_t prefix = Seq2Int(seq, n_overlap, true),
                 suffix = Seq2Int(seq.substr(seq.size() - n_overlap), n_overlap, true);
        bool is_prefix_rc(false), is_suffix_rc(false);
        if (!stranded)
        {
            uint64_t prefix_rc = GetRC(prefix, n_overlap), suffix_rc = GetRC(suffix, n_overlap);
            if (prefix_rc < prefix)
            {
                prefix = prefix_rc;
                is_prefix_rc = true;
            }
            if (suffix_rc < suffix)
            {
                suffix = suffix_rc;
                is_suffix_rc = true;
            }
        }
        {
            auto iter = hashed_mergeknot_list.insert({prefix, MergeKnot()}).first;
            std::string which_to_set = (is_prefix_rc ? "pred" : "succ");
            iter->second.AddContig(contig_serial, is_prefix_rc, which_to_set);
        }
        {
            auto iter = hashed_mergeknot_list.insert({suffix, MergeKnot()}).first;
            std::string which_to_set = (is_suffix_rc ? "succ" : "pred");
            iter->second.AddContig(contig_serial, is_suffix_rc, which_to_set);
        }
    }
}

const bool DoExtension(contigvect_t &contig_vect,
                       const fix2knot_t &hashed_mergeknot_list,
                       const CountTab &kmer_count_tab,
                       const size_t n_overlap,
                       const std::string &interv_method, const float interv_thres,
                       std::ifstream &idx_file)
{
    bool has_new_extensions(false);
    for (const auto &mk : hashed_mergeknot_list)
    {
        if (!mk.second.IsMergeable())
        {
            continue;
        }
        const uint64_t pred_serial = mk.second.GetSerial("pred"),
                       succ_serial = mk.second.GetSerial("succ");
        if (contig_vect[pred_serial].IsUsed() || contig_vect[succ_serial].IsUsed())
        {
            continue;
        }
        bool is_pred_rc = mk.second.IsRC("pred"), is_succ_rc = mk.second.IsRC("succ");
        if (interv_method != "none")
        {
            size_t left_adj_serial = contig_vect[pred_serial].GetRearKMerSerial(is_pred_rc),
                   right_adj_serial = contig_vect[succ_serial].GetHeadKMerSerial(is_succ_rc);
            if (kmer_count_tab.CalcCountDistance(left_adj_serial, right_adj_serial, interv_method, idx_file) >= interv_thres)
            {
                continue;
            }
        }
        // merge by guaranting representative k-mer having minimum p-value or input order //
        if (contig_vect[pred_serial].GetScore("origin") <= contig_vect[succ_serial].GetScore("origin")) // merge right to left
        {
            if (is_pred_rc) // prevent base contig from reverse-complement transformation, for being coherent with merging knot
            {
                contig_vect[pred_serial].LeftExtend(contig_vect[succ_serial], !is_succ_rc, n_overlap);
            }
            else
            {
                contig_vect[pred_serial].RightExtend(contig_vect[succ_serial], is_succ_rc, n_overlap);
            }
            contig_vect[succ_serial].SetUsed();
        }
        else // merge left to right
        {
            if (is_succ_rc) // prevent base contig from reverse-complement transformation, for being coherent with merging knot
            {
                contig_vect[succ_serial].RightExtend(contig_vect[pred_serial], !is_pred_rc, n_overlap);
            }
            else
            {
                contig_vect[succ_serial].LeftExtend(contig_vect[pred_serial], is_pred_rc, n_overlap);
            }
            contig_vect[pred_serial].SetUsed();
        }
        has_new_extensions = true;
    }
    return has_new_extensions;
}

void PrintContigList(const contigvect_t &contig_vect,
                     const CountTab &kmer_count_tab, const code2serial_t &code2serial,
                     const size_t k_len, const bool stranded, const size_t min_overlap,
                     const std::string &interv_method, const float interv_thres,
                     const std::string &quant_mode,
                     std::ifstream &idx_file)
{
    std::cout << "contig\tnb_merged_kmers";
    for (size_t i(0); i < kmer_count_tab.GetNbColumn(); ++i)
    {
        std::cout << "\t" << kmer_count_tab.GetColName(i);
    }
    std::cout << std::endl;
    for (const auto &elem : contig_vect)
    {
        size_t rep_uniqcode = elem.GetUniqCode(), nb_kmer_found, rep_serial = code2serial.find(rep_uniqcode)->second;
        std::string contig_seq = elem.GetSeq(), rep_kmer_seq;
        Int2Seq(rep_kmer_seq, rep_uniqcode, k_len);

        std::vector<float> sample_count;
        if (quant_mode == "rep")
        {
            kmer_count_tab.GetCountVect(sample_count, rep_serial, idx_file);
        }
        else if (quant_mode == "mean")
        {
            nb_kmer_found = EstimateAvgCount(sample_count, kmer_count_tab, code2serial, contig_seq, k_len, stranded, k_len - min_overlap,
                                             interv_method, interv_thres, idx_file);
            if (nb_kmer_found != elem.GetNbKMer()) // for debug, should not happen
            {
                throw std::domain_error("member k-mer number not coherent");
            }
        }
        std::cout << contig_seq << "\t" << elem.GetNbKMer() << "\t" << rep_kmer_seq;
        for (size_t i(1); i < kmer_count_tab.GetNbColumn(); ++i)
        {
            size_t col_serial = kmer_count_tab.GetColSerial(i);
            if (kmer_count_tab.GetColNature(i) >= 0) // count column => output according to quant_mode
            {
                std::cout << "\t" << sample_count.at(col_serial);
            }
            else if (kmer_count_tab.GetColNature(i) == -1 || kmer_count_tab.GetColNature(i) == -2) // value column => output that related with rep-k-mer
            {
                std::cout << "\t" << kmer_count_tab.GetValue(rep_serial, col_serial);
            }
        }
        std::cout << std::endl;
    }
}

int MergeMain(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string kmer_count_path, smp_info_path, interv_method("none"), quant_mode("rep"), idx_path, rep_colname;
    float interv_thres(1);
    size_t k_len(0);
    bool stranded(true);
    size_t min_overlap(0);

    ParseOptions(argc, argv, k_len, stranded, min_overlap, smp_info_path, interv_method, interv_thres, quant_mode, rep_colname, idx_path, kmer_count_path);
    PrintRunInfo(k_len, stranded, min_overlap, smp_info_path, interv_method, interv_thres, quant_mode, rep_colname, idx_path, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    contigvect_t contig_vect;          // list of contigs for extension
    code2serial_t code2serial;         // from contig's representative k-mer code to serial number in contig list
    CountTab kmer_count_tab(idx_path); // count table containing k-mer counts

    ScanCountTable(kmer_count_tab, contig_vect, code2serial, k_len, stranded, kmer_count_path, rep_colname, smp_info_path);

    std::cerr << "Count table Scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::ifstream idx_file;
    if (!idx_path.empty())
    {
        idx_file.open(idx_path);
        if (!idx_file.is_open())
        {
            throw std::domain_error("k-mer count index file " + idx_path + " was not found");
        }
    }
    for (size_t n_overlap(k_len - 1); n_overlap >= min_overlap; --n_overlap)
    {
        std::cerr << "Merging contigs with overlap " << n_overlap << std::endl;
        bool has_new_extensions(true);
        while (has_new_extensions)
        {
            std::cerr << "\tcontig list size: " << contig_vect.size() << std::endl;
            fix2knot_t hashed_mergeknot_list;
            MakeOverlapKnotDict(hashed_mergeknot_list, contig_vect, stranded, n_overlap);
            // for (const auto &elem : hashed_mergeknot_list)
            // {
            //     if (elem.second.IsMergeable())
            //     {
            //         std::string fix, contig_pred = contig_vect[elem.second.GetSerial("pred")].GetSeq(), contig_succ = contig_vect[elem.second.GetSerial("succ")].GetSeq();
            //         Int2Seq(fix, elem.first, k_len);
            //         std::cout << fix << ": " << contig_pred << " ======= " << contig_succ << std::endl;
            //     }
            // }
            has_new_extensions = DoExtension(contig_vect, hashed_mergeknot_list, kmer_count_tab, n_overlap, interv_method, interv_thres, idx_file);
            contig_vect.erase(std::remove_if(contig_vect.begin(), contig_vect.end(),
                                             [](const ContigElem &elem) { return elem.IsUsed(); }),
                              contig_vect.end());
            // PrintContigList(contig_vect, kmer_count_tab, code2serial, k_len, stranded, min_overlap, interv_method, interv_thres, quant_mode, idx_file);
        }
    }
    PrintContigList(contig_vect, kmer_count_tab, code2serial, k_len, stranded, min_overlap, interv_method, interv_thres, quant_mode, idx_file);
    if (idx_file.is_open())
    {
        idx_file.close();
    }

    std::cerr << "Contig extension finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
