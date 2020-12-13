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

void ScanCountTable(countTab_t &count_tab, TabHeader &tab_header,
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
    std::getline(kmer_count_instream, line);
    std::istringstream conv(line);
    tab_header.MakeColumnInfo(conv, "");

    std::ofstream idx_file;
    if (!idx_path.empty())
    {
        idx_file.open(idx_path);
        if (!idx_file.is_open()) // to ensure the file is opened
        {
            throw std::domain_error("error open file: " + idx_path);
        }
    }
    for (size_t iline(0); std::getline(kmer_count_instream, line); ++iline)
    {
        conv.str(line);
        count_tab.emplace_back(conv, idx_file, tab_header);
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
        if (rep_colname.empty())
        {
            contig_vect.emplace_back(std::move(seq), iline, kmer_uniqcode, iline);
        }
        else
        {
            contig_vect.emplace_back(std::move(seq), count_tab[iline].GetValue(), kmer_uniqcode, iline);
        }
        conv.clear();
    }
    if (idx_file.is_open())
    {
        idx_file.close();
    }
    count_tab.shrink_to_fit();
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
        hashed_mergeknot_list.insert({prefix, MergeKnot()}).first->second.AddContig(contig_serial, is_prefix_rc, (is_prefix_rc ? "pred" : "succ"));
        hashed_mergeknot_list.insert({suffix, MergeKnot()}).first->second.AddContig(contig_serial, is_suffix_rc, (is_suffix_rc ? "succ" : "pred"));
    }
}

const bool DoExtension(contigvect_t &contig_vect,
                       const fix2knot_t &hashed_mergeknot_list,
                       countTab_t &count_tab,
                       const size_t n_overlap,
                       const std::string &interv_method, const float interv_thres,
                       std::ifstream &idx_file, const size_t nb_counts)
{
    static std::vector<float> pred_count_vect, succ_count_vect;
    bool has_new_extensions(false);
    for (const auto &mk : hashed_mergeknot_list)
    {
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
        const TabElem &pred_feature_elem = count_tab[pred_contig.GetRearKMerSerial(is_pred_rc)],
                      &succ_feature_elem = count_tab[succ_contig.GetHeadKMerSerial(is_succ_rc)];

        if (interv_method != "none" &&
            CalcXDist(pred_feature_elem.GetCountVect(pred_count_vect, idx_file, nb_counts),
                      succ_feature_elem.GetCountVect(succ_count_vect, idx_file, nb_counts), interv_method) >= interv_thres)
        {
            continue;
        }
        // merge by guaranting representative k-mer having minimum p-value or input order //
        if (pred_contig.GetRepValue() <= succ_contig.GetRepValue()) // merge right to left
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
                     TabHeader &tab_header, countTab_t &count_tab, const code2serial_t &code2serial,
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

    std::string output_row_str;
    std::cout << "contig\tnb_merged_kmers\t" << tab_header.MakeOutputHeaderStr(output_row_str) << std::endl;
    for (size_t i(0); i < tab_header.GetNbCol(); ++i)
    {
        std::cout << "\t" << tab_header.GetColNameAt(i);
    }
    std::cout << std::endl;

    std::string rep_kmer_seq;
    size_t rep_uniqcode, rep_serial, nb_count = tab_header.GetNbCount();
    if (quant_mode == "rep")
    {
        for (const auto &elem : contig_vect)
        {
            rep_uniqcode = elem.GetRepUniqcode();
            rep_serial = code2serial.find(rep_uniqcode)->second;
            Int2Seq(rep_kmer_seq, rep_uniqcode, k_len);
            std::cout << elem.GetSeq() << "\t" << elem.GetNbKMer() << "\t" << rep_kmer_seq << "\t"
                      << count_tab[rep_serial].MakeOutputRowStr(output_row_str, idx_file, nb_count) << std::endl;
        }
    }
    else if (quant_mode == "mean" || quant_mode == "median")
    {
        std::vector<float> count_vect;
        std::vector<std::vector<float>> count_vect_vect;
        for (const auto &elem : contig_vect)
        {
            rep_uniqcode = elem.GetRepUniqcode();
            rep_serial = code2serial.find(rep_uniqcode)->second;
            Int2Seq(rep_kmer_seq, rep_uniqcode, k_len);
            count_vect_vect.resize(nb_count, std::vector<float>(elem.GetNbKMer()));
            for (size_t rs : elem.GetMemKMerSerialVect())
            {
                count_tab[rs].GetCountVect(count_vect, idx_file, nb_count);
                for (size_t i(0); i < nb_count; ++i)
                {
                    count_vect_vect[i].emplace_back(count_vect[i]);
                }
            }
            std::cout << elem.GetSeq() << "\t" << elem.GetNbKMer() << "\t" << rep_kmer_seq << "\t"
                      << count_tab[rep_serial].MakeOutputRowStr(output_row_str, CalcVectsRes(count_vect, count_vect_vect, quant_mode), idx_file)
                      << std::endl;
        }
    }
    else
    {
        throw std::domain_error("unknown quant mode: " + quant_mode);
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
    countTab_t count_tab;                // feature count table

    ScanCountTable(count_tab, tab_header, contig_vect, code2serial, k_len, stranded, kmer_count_path, rep_colname, idx_path);

    std::cerr << "Count table scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
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
            has_new_extensions = DoExtension(contig_vect, hashed_mergeknot_list, count_tab, n_overlap, interv_method, interv_thres, idx_file, tab_header.GetNbCount());
            contig_vect.erase(std::remove_if(contig_vect.begin(), contig_vect.end(),
                                             [](const ContigElem &elem) { return elem.IsUsed(); }),
                              contig_vect.end());
            contig_vect.shrink_to_fit();
            // PrintContigList(contig_vect, kmer_count_tab, code2serial, k_len, stranded, min_overlap, interv_method, interv_thres, quant_mode, idx_file);
        }
    }

    std::cerr << "Contig extension finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    PrintContigList(contig_vect, tab_header, count_tab, code2serial, k_len, quant_mode, idx_file, out_path);
    if (idx_file.is_open())
    {
        idx_file.close();
    }

    std::cerr << "Contig print finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
