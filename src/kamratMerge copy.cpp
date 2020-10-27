#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include <stdexcept>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "utils/seq_coding.hpp"
#include "utils/vec_operation.hpp"
#include "data_struct/feature_tab_elem.hpp"
#include "data_struct/feature_tab_header.hpp"
#include "data_struct/contig_elem.hpp"
#include "data_struct/merge_knot.hpp"
#include "run_info_parser/merge.hpp"
#include "run_info_parser/utils.hpp"

void ScanCountTable(featuretab_t &feature_tab, FeatureTabHeader &feature_tab_header,
                    contigvect_t &contig_vect, code2serial_t &code2serial,
                    const size_t k_len, const bool stranded,
                    const std::string &kmer_count_path,
                    const std::string &rep_colname,
                    const std::string &idx_path)
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

    std::string line, seq;
    //----- Dealing with Header Line for Constructing ColumnInfo Object -----//
    std::getline(kmer_count_instream, line);
    size_t rep_colpos = feature_tab_header.MakeColumnInfo(line, rep_colname),
           nb_value = feature_tab_header.GetNbValue(), nb_count = feature_tab_header.GetNbCount(), nb_str = feature_tab_header.GetNbStr();
    //----- Dealing with Following k-mer Count Lines -----//
    std::ofstream count_idx_file;
    if (!idx_path.empty())
    {
        count_idx_file.open(idx_path);
        if (!count_idx_file.is_open()) // to ensure the file is opened
        {
            throw std::domain_error("error open file: " + idx_path);
        }
    }
    std::istringstream conv;
    for (size_t iline(0); std::getline(kmer_count_instream, line); ++iline)
    {
        conv.str(line);
        feature_tab.emplace_back(nb_value, nb_count, nb_str, conv, count_idx_file, feature_tab_header.GetColNatureVect());
        float rep_val = (rep_colpos == 0 ? iline : feature_tab.back().GetValueAt(feature_tab_header.GetColSerialVect().at(rep_colpos)));
        seq = std::move(line.substr(0, line.find_first_of(" \t"))); // first column as feature (string)
        if (seq.size() != k_len)
        {
            throw std::domain_error("the given k-len parameter not coherent with input k-mer: " + seq);
        }
        uint64_t kmer_uniqcode = Seq2Int(seq, seq.size(), stranded);
        if (!code2serial.insert({kmer_uniqcode, iline}).second)
        {
            throw std::domain_error("duplicate input (newly inserted k-mer already exists in k-mer hash list)");
        }
        contig_vect.emplace_back(std::move(seq), rep_val, kmer_uniqcode, iline);
    }
    if (count_idx_file.is_open())
    {
        count_idx_file.close();
    }
    feature_tab.shrink_to_fit();
    contig_vect.shrink_to_fit();
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
                       featuretab_t &feature_count_tab,
                       const size_t n_overlap,
                       const std::string &interv_method, const float interv_thres,
                       std::ifstream &idx_file, const size_t nb_count)
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
        if (idx_file.is_open())
        {
            feature_count_tab[contig_vect[pred_serial].GetRearKMerSerial(is_pred_rc)].RestoreCountVect(idx_file, nb_count);
            feature_count_tab[contig_vect[succ_serial].GetHeadKMerSerial(is_succ_rc)].RestoreCountVect(idx_file, nb_count);
        }
        if (interv_method != "none" &&
            CalcXDist(feature_count_tab[contig_vect[pred_serial].GetRearKMerSerial(is_pred_rc)].GetCountVect(),
                      feature_count_tab[contig_vect[succ_serial].GetHeadKMerSerial(is_succ_rc)].GetCountVect(),
                      interv_method) >= interv_thres)
        {
            feature_count_tab[contig_vect[pred_serial].GetRearKMerSerial(is_pred_rc)].ClearCountVect();
            feature_count_tab[contig_vect[succ_serial].GetHeadKMerSerial(is_succ_rc)].ClearCountVect();
            continue;
        }
        feature_count_tab[contig_vect[pred_serial].GetRearKMerSerial(is_pred_rc)].ClearCountVect();
        feature_count_tab[contig_vect[succ_serial].GetHeadKMerSerial(is_succ_rc)].ClearCountVect();
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
                     const std::string &interv_method, const float interv_thres,
                     const std::string &quant_mode,
                     std::ifstream &idx_file,
                     const std::string &out_path)
{
    std::ofstream out_file;
    if (!out_path.empty())
    {
        out_file.open(out_path);
        if (!out_file.is_open())
        {
            throw std::domain_error("cannot open file: " + out_path);
        }
    }
    auto backup_buf = std::cout.rdbuf();
    if (!out_path.empty()) // output to file if a path is given, to screen if not
    {
        std::cout.rdbuf(out_file.rdbuf());
    }
    const auto k_len = kmer_count_tab.GetKLen();
    std::cout << "contig\tnb_merged_kmers";
    for (size_t i(0); i < kmer_count_tab.GetNbColumn(); ++i)
    {
        std::cout << "\t" << kmer_count_tab.GetColName(i);
    }
    std::cout << std::endl;

    std::string contig_seq, rep_kmer_seq;
    std::vector<float> sample_count;
    for (const auto &elem : contig_vect)
    {
        size_t rep_uniqcode = elem.GetUniqCode(), rep_serial = code2serial.find(rep_uniqcode)->second;
        contig_seq = elem.GetSeq();
        Int2Seq(rep_kmer_seq, rep_uniqcode, k_len);

        if (quant_mode == "rep")
        {
            kmer_count_tab.GetCountVect(sample_count, rep_serial, idx_file);
        }
        else if (quant_mode == "mean")
        {
            kmer_count_tab.EstimateMeanCountVect(sample_count, elem.GetMemKMerSerialVect(), idx_file);
        }
        else
        {
            throw std::domain_error("unknown quant mode: " + quant_mode);
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
        rep_kmer_seq.clear();
        sample_count.clear();
    }
    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
}

int MergeMain(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string kmer_count_path, smp_info_path, interv_method("none"), quant_mode("rep"), idx_path, out_path, rep_colname;
    float interv_thres(1);
    size_t k_len(0);
    bool stranded(true);
    size_t min_overlap(0);

    ParseOptions(argc, argv, k_len, stranded, min_overlap, smp_info_path, interv_method, interv_thres, quant_mode, rep_colname, idx_path, out_path, kmer_count_path);
    PrintRunInfo(k_len, stranded, min_overlap, smp_info_path, interv_method, interv_thres, quant_mode, rep_colname, idx_path, out_path, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    contigvect_t contig_vect;                           // list of contigs for extension
    code2serial_t code2serial;                          // from contig's representative k-mer code to serial number in contig list
    CountTab kmer_count_tab(k_len, stranded, idx_path); // count table containing k-mer counts

    ScanCountTable(kmer_count_tab, contig_vect, code2serial, kmer_count_path, rep_colname, smp_info_path);

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
    for (size_t n_overlap = k_len - 1; n_overlap >= min_overlap; --n_overlap)
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
            contig_vect.shrink_to_fit();
            // PrintContigList(contig_vect, kmer_count_tab, code2serial, k_len, stranded, min_overlap, interv_method, interv_thres, quant_mode, idx_file);
        }
    }

    std::cerr << "Contig extension finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    PrintContigList(contig_vect, kmer_count_tab, code2serial, interv_method, interv_thres, quant_mode, idx_file, out_path);
    if (idx_file.is_open())
    {
        idx_file.close();
    }

    std::cerr << "Contig print finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
