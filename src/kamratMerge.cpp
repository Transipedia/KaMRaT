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
#include "data_struct/tab_elem.hpp"
#include "data_struct/tab_header.hpp"
#include "data_struct/contig_elem.hpp"
#include "data_struct/merge_knot.hpp"
#include "run_info_parser/merge.hpp"
#include "run_info_parser/utils.hpp"

void ScanCountTable(featuretab_t &feature_tab, TabHeader &tab_header,
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
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    {
        size_t pos = kmer_count_path.find_last_of(".");
        if (pos != std::string::npos && kmer_count_path.substr(pos + 1) == "gz")
        {
            inbuf.push(boost::iostreams::gzip_decompressor());
        }
    }
    inbuf.push(kmer_count_file);
    std::istream kmer_count_instream(&inbuf);

    std::string line, seq;
    //----- Dealing with Header Line for Constructing ColumnInfo Object -----//
    std::getline(kmer_count_instream, line);
    tab_header.MakeColumnInfo(line, rep_colname);
    //----- Dealing with Following k-mer Count Lines -----//
    std::ofstream idx_file;
    if (!idx_path.empty())
    {
        idx_file.open(idx_path);
        if (!idx_file.is_open()) // to ensure the file is opened
        {
            throw std::domain_error("error open file: " + idx_path);
        }
    }
    std::istringstream conv;
    float rep_val;
    for (size_t iline(0); std::getline(kmer_count_instream, line); ++iline)
    {
        conv.str(line);
        feature_tab.emplace_back(conv, idx_file, rep_val, tab_header);
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
        conv.clear();
    }
    if (idx_file.is_open())
    {
        idx_file.close();
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
    static std::string seq;
    for (size_t contig_serial(0); contig_serial < contig_vect.size(); ++contig_serial)
    {
        seq = contig_vect[contig_serial].GetSeq();
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
            iter->second.AddContig(contig_serial, is_prefix_rc, (is_prefix_rc ? "pred" : "succ"));
        }
        {
            auto iter = hashed_mergeknot_list.insert({suffix, MergeKnot()}).first;
            iter->second.AddContig(contig_serial, is_suffix_rc, (is_suffix_rc ? "succ" : "pred"));
        }
    }
}

const bool DoExtension(contigvect_t &contig_vect,
                       const fix2knot_t &hashed_mergeknot_list,
                       featuretab_t &feature_tab,
                       const size_t n_overlap,
                       const std::string &interv_method, const float interv_thres,
                       std::ifstream &idx_file, const TabHeader &tab_header)
{
    static std::vector<float> pred_count_vect, succ_count_vect;
    const size_t nb_counts = tab_header.GetNbCount();
    bool has_new_extensions(false);
    for (const auto &mk : hashed_mergeknot_list)
    {
        pred_count_vect.clear();
        succ_count_vect.clear();
        
        if (!mk.second.IsMergeable())
        {
            continue;
        }
        ContigElem &pred_contig = contig_vect[mk.second.GetSerial("pred")],
                   &succ_contig = contig_vect[mk.second.GetSerial("succ")];
        if (pred_contig.IsUsed() || succ_contig.IsUsed())
        {
            continue;
        }
        bool is_pred_rc = mk.second.IsRC("pred"), is_succ_rc = mk.second.IsRC("succ");
        const TabElem &pred_feature_elem = feature_tab[pred_contig.GetRearKMerSerial(is_pred_rc)],
                      &succ_feature_elem = feature_tab[succ_contig.GetHeadKMerSerial(is_succ_rc)];

        if (interv_method != "none" &&
            CalcXDist(pred_feature_elem.GetCountVect(pred_count_vect, idx_file, nb_counts),
                      succ_feature_elem.GetCountVect(succ_count_vect, idx_file, nb_counts), interv_method) >= interv_thres)
        {
            continue;
        }
        // merge by guaranting representative k-mer having minimum p-value or input order //
        if (pred_contig.GetScore("origin") <= succ_contig.GetScore("origin")) // merge right to left
        {
            if (is_pred_rc) // prevent base contig from reverse-complement transformation, for being coherent with merging knot
            {
                pred_contig.LeftExtend(succ_contig, !is_succ_rc, n_overlap);
            }
            else
            {
                pred_contig.RightExtend(succ_contig, is_succ_rc, n_overlap);
            }
            succ_contig.SetUsed();
        }
        else // merge left to right
        {
            if (is_succ_rc) // prevent base contig from reverse-complement transformation, for being coherent with merging knot
            {
                succ_contig.RightExtend(pred_contig, !is_pred_rc, n_overlap);
            }
            else
            {
                succ_contig.LeftExtend(pred_contig, is_pred_rc, n_overlap);
            }
            pred_contig.SetUsed();
        }
        has_new_extensions = true;
    }
    return has_new_extensions;
}

void PrintContigList(const contigvect_t &contig_vect,
                     TabHeader &tab_header, featuretab_t &feature_tab, const code2serial_t &code2serial,
                     const size_t k_len, const std::string &quant_mode,
                     std::ifstream &idx_file, const std::string &out_path)
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
    std::cout << "contig\tnb_merged_kmers";
    for (size_t i(0); i < tab_header.GetNbCol(); ++i)
    {
        std::cout << "\t" << tab_header.GetColNameAt(i);
    }
    std::cout << std::endl;

    std::string rep_kmer_seq;
    size_t rep_uniqcode, rep_serial, nb_count = tab_header.GetNbCount(), nb_value = tab_header.GetNbValue();
    std::vector<float> count_vect, value_vect, count_vect_x, value_vect_x;
    for (const auto &elem : contig_vect)
    {
        rep_uniqcode = elem.GetUniqCode();
        rep_serial = code2serial.find(rep_uniqcode)->second;
        Int2Seq(rep_kmer_seq, rep_uniqcode, k_len);

        if (quant_mode == "rep")
        {
            feature_tab[rep_serial].GetVectsAndClear(count_vect, value_vect, idx_file, nb_count, nb_value);
        }
        else if (quant_mode == "mean")
        {
            count_vect.assign(nb_count, 0);
            for (size_t rs : elem.GetMemKMerSerialVect())
            {
                feature_tab[rs].GetVectsAndClear(count_vect_x, value_vect_x, idx_file, nb_count, nb_value);
                for (size_t col_serial(0); col_serial < nb_count; ++col_serial)
                {
                    count_vect[col_serial] += count_vect_x[col_serial];
                }
                if (rs == rep_serial)
                {
                    value_vect = std::move(value_vect_x);
                }
            }
            for (size_t col_serial(0); col_serial < nb_count; ++col_serial)
            {
                count_vect[col_serial] /= elem.GetMemKMerSerialVect().size();
            }
        }
        else
        {
            throw std::domain_error("unknown quant mode: " + quant_mode);
        }
        std::cout << elem.GetSeq() << "\t" << elem.GetNbKMer() << "\t" << rep_kmer_seq;
        for (size_t i(1); i < tab_header.GetNbCol(); ++i)
        {
            size_t col_serial = tab_header.GetColSerialAt(i);
            if (tab_header.GetColNatureAt(i) >= 0) // count column => output according to quant_mode
            {
                std::cout << "\t" << count_vect[col_serial];
            }
            else if (tab_header.GetColNatureAt(i) == -1 || tab_header.GetColNatureAt(i) == -2) // value column => output that related with rep-k-mer
            {
                std::cout << "\t" << value_vect[col_serial];
            }
        }
        std::cout << std::endl;
        rep_kmer_seq.clear();
        count_vect_x.clear();
        value_vect_x.clear();
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

    contigvect_t contig_vect;            // list of contigs for extension
    code2serial_t code2serial;           // from contig's representative k-mer code to serial number in contig list
    TabHeader tab_header(smp_info_path); // the header of feature table
    featuretab_t feature_tab;            // feature count table

    ScanCountTable(feature_tab, tab_header, contig_vect, code2serial, k_len, stranded, kmer_count_path, rep_colname, idx_path);

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
            has_new_extensions = DoExtension(contig_vect, hashed_mergeknot_list, feature_tab, n_overlap, interv_method, interv_thres, idx_file, tab_header);
            contig_vect.erase(std::remove_if(contig_vect.begin(), contig_vect.end(),
                                             [](const ContigElem &elem) { return elem.IsUsed(); }),
                              contig_vect.end());
            contig_vect.shrink_to_fit();
            // PrintContigList(contig_vect, kmer_count_tab, code2serial, k_len, stranded, min_overlap, interv_method, interv_thres, quant_mode, idx_file);
        }
    }

    std::cerr << "Contig extension finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    PrintContigList(contig_vect, tab_header, feature_tab, code2serial, k_len, quant_mode, idx_file, out_path);
    if (idx_file.is_open())
    {
        idx_file.close();
    }

    std::cerr << "Contig print finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
