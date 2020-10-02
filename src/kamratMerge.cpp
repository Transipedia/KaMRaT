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

#include "utils/statistics.hpp"
#include "utils/seq_coding.hpp"
#include "data_struct/count_tab_by_fields.hpp"
#include "data_struct/contig_elem.hpp"
#include "data_struct/merge_knot.hpp"
#include "run_info_parser/merge.hpp"
#include "run_info_parser/utils.hpp"

void ScanCountTable(CountTabByFields &kmer_count_tab,
                    code2contig_t &hashed_contig_list,
                    const std::string &kmer_count_path,
                    const std::string &score_colname,
                    const std::string &sample_info_path,
                    const std::string &count_index_path,
                    const bool stranded)
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
    kmer_count_tab.MakeSmpCond(sample_info_path);
    kmer_count_tab.MakeColumnInfo(line, score_colname);
    //----- Dealing with Following k-mer Count Lines -----//
    std::ofstream count_index_file(count_index_path);
    for (size_t kmer_serial(0); std::getline(kmer_count_instream, line); ++kmer_serial)
    {
        // std::cout << line << std::endl;
        float rep_val;
        bool has_rep_val;
        if (kmer_count_tab.GetMode() == "inMem")
        {
            has_rep_val = kmer_count_tab.AddCountInMem(rep_val, line);
        }
        else if (kmer_count_tab.GetMode() == "onDsk")
        {
            has_rep_val = kmer_count_tab.AddIndexOnDsk(rep_val, line, count_index_file);
        }
        if (!has_rep_val) // if no representative k-mer column
        {
            rep_val = kmer_serial;
        }
        std::istringstream conv(line);
        std::string seq;
        conv >> seq; // first column as feature (string)
        if (!hashed_contig_list.insert({Seq2Int(seq, seq.size(), stranded), ContigElem(seq, rep_val, kmer_serial)}).second)
        {
            throw std::domain_error("duplicate input (newly inserted k-mer already exists in k-mer hash list)");
        }
    }
    count_index_file.close();
    kmer_count_file.close();
}

void MakeOverlapKnotDict(fix2knot_t &hashed_mergeknot_list,
                         const code2contig_t &hashed_contig_list,
                         const bool stranded,
                         const size_t n_overlap)
{
    for (const auto &elem : hashed_contig_list)
    {
        std::string seq = elem.second.GetSeq();
        if (seq.size() <= n_overlap)
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
            iter->second.AddContig(elem.first, is_prefix_rc, which_to_set);
        }
        {
            auto iter = hashed_mergeknot_list.insert({suffix, MergeKnot()}).first;
            std::string which_to_set = (is_suffix_rc ? "succ" : "pred");
            iter->second.AddContig(elem.first, is_suffix_rc, which_to_set);
        }
    }
}

void PrintContigList(const code2contig_t &hashed_contig_list,
                     const CountTabByFields &kmer_count_tab,
                     const size_t k_len,
                     const std::string &quant_mode,
                     std::ifstream &index_file)
{
    std::cout << "contig\tnb_merged_kmers";
    for (size_t i(0); i < kmer_count_tab.GetNbColumn(); ++i)
    {
        std::cout << "\t" << kmer_count_tab.GetColName(i);
    }
    std::cout << std::endl;
    for (const auto &elem : hashed_contig_list)
    {
        std::string rep_kmer_seq;
        size_t rep_serial = elem.second.GetRepSerial();
        Int2Seq(rep_kmer_seq, elem.first, k_len);
        std::cout << elem.second.GetSeq() << "\t" << elem.second.GetNbMemberKMer() << "\t" << rep_kmer_seq;
        std::vector<float> sample_count;
        if (quant_mode == "rep" && kmer_count_tab.GetMode() == "inMem")
        {
            kmer_count_tab.GetCountInMem(sample_count, rep_serial);
        }
        else if (quant_mode == "rep" && kmer_count_tab.GetMode() == "onDsk")
        {
            kmer_count_tab.GetCountOnDsk(sample_count, rep_serial, index_file);
        }
        else if (quant_mode == "mean" && kmer_count_tab.GetMode() == "inMem")
        {
            kmer_count_tab.GetAvgCountInMem(sample_count, elem.second.GetKMerSerialSet());
        }
        else if (quant_mode == "mean" && kmer_count_tab.GetMode() == "onDsk")
        {
            kmer_count_tab.GetAvgCountOnDsk(sample_count, elem.second.GetKMerSerialSet(), index_file);
        }
        for (size_t i(1); i < kmer_count_tab.GetNbColumn(); ++i)
        {
            size_t serial = kmer_count_tab.GetColSerial(i);
            if (kmer_count_tab.GetColNature(i) >= 0) // count column => output according to quant_mode
            {
                std::cout << "\t" << sample_count.at(serial);
            }
            else if (kmer_count_tab.GetColNature(i) == -1 || kmer_count_tab.GetColNature(i) == -2) // value column => output that related with rep-k-mer
            {
                std::cout << "\t" << kmer_count_tab.GetValue(rep_serial, serial);
            }
        }
        std::cout << std::endl;
    }
}

const bool DoExtension(code2contig_t &hashed_contig_list,
                       const fix2knot_t &hashed_mergeknot_list,
                       const CountTabByFields &kmer_count_tab,
                       const size_t n_overlap,
                       const std::string &interv_method,
                       const float interv_thres,
                       std::ifstream &index_file)
{
    bool has_new_extensions(false);
    for (const auto &mk : hashed_mergeknot_list)
    {
        if (!mk.second.IsMergeable())
        {
            continue;
        }
        const uint64_t pred_code(mk.second.GetPredCode()), succ_code(mk.second.GetSuccCode());
        const auto pred_iter = hashed_contig_list.find(pred_code), succ_iter = hashed_contig_list.find(succ_code);
        if (pred_iter == hashed_contig_list.cend() || succ_iter == hashed_contig_list.cend())
        {
            continue;
        }
        if (interv_method != "none")
        {
            size_t left_serial = (mk.second.IsPredRC() ? pred_iter->second.GetHeadSerial() : pred_iter->second.GetRearSerial()),
                   right_serial = (mk.second.IsSuccRC() ? succ_iter->second.GetRearSerial() : succ_iter->second.GetHeadSerial());
            std::vector<float> pred_counts, succ_counts;
            if (kmer_count_tab.GetMode() == "onDsk")
            {
                kmer_count_tab.GetCountOnDsk(pred_counts, left_serial, index_file);
                kmer_count_tab.GetCountOnDsk(succ_counts, right_serial, index_file);
            }
            else
            {
                kmer_count_tab.GetCountInMem(pred_counts, left_serial);
                kmer_count_tab.GetCountInMem(succ_counts, right_serial);
            }
            if (interv_method == "pearson" && CalcPearsonCorrelation(pred_counts, succ_counts) < interv_thres)
            {
                continue;
            }
            else if (interv_method == "spearman" && CalcSpearmanCorrelation(pred_counts, succ_counts) < interv_thres)
            {
                continue;
            }
            else if (interv_method == "mac" && CalcMeanAbsoluteContrast(pred_counts, succ_counts) > interv_thres)
            {
                continue;
            }
        }
        if (mk.second.IsPredRC())
        {
            pred_iter->second.SelfReverseComplement();
        }
        if (mk.second.IsSuccRC())
        {
            succ_iter->second.SelfReverseComplement();
        }
        if (pred_iter->second.GetScore("origin") <= succ_iter->second.GetScore("origin"))
        {
            pred_iter->second.RightMerge(succ_iter->second, n_overlap);
            hashed_contig_list.erase(succ_iter);
            if (mk.second.IsPredRC())
            {
                pred_iter->second.SelfReverseComplement();
            }
        }
        else // pred_iter->second.GetRepValue() > succ_iter->second.GetRepValue()
        {
            succ_iter->second.LeftMerge(pred_iter->second, n_overlap);
            hashed_contig_list.erase(pred_iter);
            if (mk.second.IsSuccRC())
            {
                succ_iter->second.SelfReverseComplement();
            }
        }
        has_new_extensions = true;
    }
    return has_new_extensions;
}

int MergeMain(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string kmer_count_path, sample_info_path, interv_method("none"), quant_mode("rep"), tmp_dir("./"), rep_colname;
    float interv_thres;
    size_t k_len(0);
    bool stranded(true), disk_mode(false);
    size_t min_overlap(0);

    ParseOptions(argc, argv, k_len, stranded, min_overlap, sample_info_path, interv_method, interv_thres, quant_mode, rep_colname, disk_mode, tmp_dir, kmer_count_path);
    PrintRunInfo(k_len, stranded, min_overlap, sample_info_path, interv_method, interv_thres, quant_mode, rep_colname, disk_mode, tmp_dir, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    code2contig_t hashed_contig_list;
    CountTabByFields kmer_count_tab((disk_mode ? "onDsk" : "inMem"));
    std::string count_index_path = tmp_dir + "counts.idx";
    ScanCountTable(kmer_count_tab, hashed_contig_list, kmer_count_path, rep_colname, sample_info_path, count_index_path, stranded);

    std::cerr << "Count table Scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::ifstream index_file(count_index_path);
    if (!index_file.is_open())
    {
        throw std::domain_error("k-mer count index file " + count_index_path + " was not found");
    }
    for (size_t n_overlap(k_len - 1); n_overlap >= min_overlap; --n_overlap)
    {
        std::cerr << "Merging contigs with overlap " << n_overlap << std::endl;
        bool has_new_extensions(true);
        while (has_new_extensions)
        {
            std::cerr << "\tcontig list size: " << hashed_contig_list.size() << std::endl;
            fix2knot_t hashed_mergeknot_list;
            MakeOverlapKnotDict(hashed_mergeknot_list, hashed_contig_list, stranded, n_overlap);
            has_new_extensions = DoExtension(hashed_contig_list, hashed_mergeknot_list, kmer_count_tab, n_overlap, interv_method, interv_thres, index_file);
        }
    }
    PrintContigList(hashed_contig_list, kmer_count_tab, k_len, quant_mode, index_file);
    index_file.close();

    std::cerr << "Contig extension finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
