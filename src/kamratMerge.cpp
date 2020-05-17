#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>

#include "utils/utils.hpp"
#include "utils/statistics.hpp"
#include "utils/column_info.hpp"
#include "utils/seq_coding.hpp"
#include "merge/contig_info.hpp"
#include "merge/parse_opt_print_info.hpp"
#include "merge/utils.hpp"
#include "merge/overlap_knot.hpp"
#include "merge/dna_sequence.hpp"

void ScanCountTable(std::vector<ContigInfo> &contig_list,
                    ColumnInfo &column_info,
                    const std::string &kmer_count_path,
                    const std::string &rep_value_cname,
                    const std::string &sample_info_path,
                    const std::string &count_index_path,
                    const bool stranded,
                    const bool disk_mode)
{
    std::ifstream kmer_count_file(kmer_count_path);
    ExitIf(!kmer_count_file.is_open(), "ERROR: k-mer count file " + kmer_count_path + " was not found.");
    std::string line;
    // dealing with the header line for parsing column names //
    std::getline(kmer_count_file, line);
    column_info.MakeColumnInfo(line, sample_info_path, rep_value_cname);
    // dealing with the count lines //
    std::ofstream count_index_file(count_index_path);
    std::unordered_set<uint64_t> kmer_set;
    for (size_t i_line(1); std::getline(kmer_count_file, line); ++i_line)
    {
        std::istringstream conv(line);
        std::string term, seq;
        std::vector<double> sample_counts;
        float rep_value;
        for (size_t i(0); conv >> term; ++i)
        {
            char col_nat = column_info.GetColumnNature(i);
            if (col_nat == 'f')
            {
                seq = term;
            }
            else if (col_nat == 's')
            {
                sample_counts.push_back(std::stof(term));
            }
            else if (col_nat == 'r')
            {
                rep_value = std::stof(term);
            }
        }
        if (rep_value_cname.empty())
        {
            rep_value = i_line;
        }
        ExitIf(!kmer_set.insert(Seq2Int(seq, seq.size(), stranded)).second, "ERROR: duplicated k-mer " + contig_list.back().GetSequence());
        if (disk_mode)
        {
            contig_list.emplace_back(seq, rep_value, WriteCountToIndex(count_index_file, sample_counts));
        }
        else
        {
            contig_list.emplace_back(seq, rep_value, sample_counts);
        }
    }
    count_index_file.close();
    kmer_count_file.close();
}

void MakeOverlapKnotDict(std::map<uint64_t, OverlapKnot> &overlap_knot_dict,
                         const std::vector<ContigInfo> &contig_list,
                         const bool stranded,
                         const size_t n_overlap)
{
    size_t nb_contigs(contig_list.size() - 1); // "- 1" for ignoring the leading placeholder contig
    for (size_t contig_serial(1); contig_serial <= nb_contigs; ++contig_serial)
    {
        ExitIf(contig_list.at(contig_serial).IsUsed(), "ERROR: used contig not deleted.");
        std::string seq = contig_list.at(contig_serial).GetSequence();
        ExitIf(seq.size() <= n_overlap, "ERROR: k-mer size less than k-length");
        uint64_t prefix(Seq2Int(seq, n_overlap, true)), suffix(Seq2Int(seq.substr(seq.size() - n_overlap), n_overlap, true));
        bool is_prefix_rc(false), is_suffix_rc(false);
        if (!stranded)
        {
            uint64_t prefix_rc(GetRC(prefix, n_overlap)), suffix_rc(GetRC(suffix, n_overlap));
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
            auto iter = overlap_knot_dict.insert({prefix, OverlapKnot()}).first;
            is_prefix_rc ? (iter->second.SetPredtig(contig_serial, true)) : (iter->second.SetSucctig(contig_serial, false));
        }
        {
            auto iter = overlap_knot_dict.insert({suffix, OverlapKnot()}).first;
            is_suffix_rc ? (iter->second.SetSucctig(contig_serial, true)) : (iter->second.SetPredtig(contig_serial, false));
        }
    }
}

void EstimateContigCount(std::vector<double> &counts,
                         const ContigInfo &contig,
                         const size_t nb_sample,
                         std::ifstream &index_file)
{
    // set count size //
    counts.resize(nb_sample, 0);
    size_t nb_count_pos = contig.GetCountPos().size();
    ExitIf(nb_count_pos == 0, "ERROR: invalid count source number.");
    for (size_t count_pos : contig.GetCountPos())
    {
        AddCountFromIndex(index_file, counts, count_pos, nb_sample);
    }
    // if (seq == "TTATGGAGCAGCCCCCTTTGAAGATCTCCAAGTAGACTTCAC")
    // {
    //     std::cerr << kmer_seq;
    //     for (size_t i(0); i < nb_sample; ++i)
    //     {
    //         std::cerr << "\t" << sum.at(i);
    //     }
    //     std::cerr << std::endl;
    // }
    // turn sum to average //
    for (size_t i = 0; i < nb_sample; ++i)
    {
        counts.at(i) /= nb_count_pos;
    }
    // if (seq == "TTATGGAGCAGCCCCCTTTGAAGATCTCCAAGTAGACTTCAC")
    //     std::cerr << sum.at(0) << "/" << nb_kmers << "=" << counts.at(0) << "\t" << std::endl;
}

void PrintContigList(std::ostream &out_s,
                     const std::vector<ContigInfo> &contig_list,
                     const size_t nb_sample,
                     const bool disk_mode,
                     std::ifstream &index_file)
{
    size_t nb_contigs(contig_list.size());
    for (size_t contig_serial(1); contig_serial < nb_contigs; ++contig_serial)
    {
        ContigInfo contig(contig_list.at(contig_serial));
        out_s << contig.GetSequence() << "\t" << contig.GetNbKmer() << "\t" << contig.GetTag();
        std::vector<double> counts;
        if (disk_mode)
        {
            EstimateContigCount(counts, contig, nb_sample, index_file);
        }
        else
        {
            counts = contig_list.at(contig_serial).GetCounts();
        }
        for (double count_x : counts)
        {
            out_s << std::fixed << std::setprecision(2) << "\t" << (count_x + 0.005);
        }
        out_s << std::endl;
    }
}

const bool DoExtension(std::vector<ContigInfo> &contig_list,
                       const std::map<uint64_t, OverlapKnot> &overlap_knot_dict,
                       const size_t n_overlap,
                       const std::string &interv_method,
                       const double interv_thres,
                       const size_t nb_sample,
                       const std::string &quant_mode,
                       const bool disk_mode,
                       const bool comp_adj,
                       std::ifstream &index_file)
{
    bool have_new_extensions(false);
    // std::cout << "===========" << n_overlap << "===========" << std::endl;
    for (const auto &ov_kt : overlap_knot_dict)
    {
        if (!ov_kt.second.IsSingleKnot())
        {
            continue;
        }
        const size_t pred_serial(ov_kt.second.GetPredtigSerial()), succ_serial(ov_kt.second.GetSucctigSerial());
        if (contig_list.at(pred_serial).IsUsed() || contig_list.at(succ_serial).IsUsed())
        {
            continue;
        }
        if (interv_method != "none")
        {
            std::vector<double> pred_counts, succ_counts;
            if (disk_mode)
            {
                if (comp_adj)
                {
                    ov_kt.second.IsPredtigRC() ? LoadCountFromIndex(index_file, pred_counts, contig_list.at(pred_serial).GetHeadPos(), nb_sample)
                                               : LoadCountFromIndex(index_file, pred_counts, contig_list.at(pred_serial).GetRearPos(), nb_sample);
                    ov_kt.second.IsSucctigRC() ? LoadCountFromIndex(index_file, succ_counts, contig_list.at(succ_serial).GetRearPos(), nb_sample)
                                               : LoadCountFromIndex(index_file, succ_counts, contig_list.at(succ_serial).GetHeadPos(), nb_sample);
                }
                else
                {
                    EstimateContigCount(pred_counts, contig_list.at(pred_serial), nb_sample, index_file);
                    EstimateContigCount(succ_counts, contig_list.at(succ_serial), nb_sample, index_file);
                }
            }
            else
            {
                if (comp_adj)
                {
                    pred_counts = ov_kt.second.IsPredtigRC() ? contig_list.at(pred_serial).GetHeadCounts() : contig_list.at(pred_serial).GetRearCounts();
                    succ_counts = ov_kt.second.IsSucctigRC() ? contig_list.at(succ_serial).GetRearCounts() : contig_list.at(succ_serial).GetHeadCounts();
                }
                else
                {
                    pred_counts = contig_list.at(pred_serial).GetCounts();
                    succ_counts = contig_list.at(succ_serial).GetCounts();
                }
            }
            // {
            //     std::string fix_seq;
            //     Int2Seq(fix_seq, ov_kt.first, n_overlap);
            //     if (fix_seq == "ATCAGAATAGCCACATTTA")
            //     {
            //         std::cerr << (ov_kt.second.IsSingleKnot() ? "ture" : "false") << '\t';
            //         std::string left_contig = contig_list.at(ov_kt.second.GetPredtigSerial()).GetSequence();
            //         if (ov_kt.second.IsPredtigRC())
            //         {
            //             ToReverseComplement(left_contig);
            //         }
            //         std::string right_contig = contig_list.at(ov_kt.second.GetSucctigSerial()).GetSequence();
            //         if (ov_kt.second.IsSucctigRC())
            //         {
            //             ToReverseComplement(right_contig);
            //         }
            //         std::cerr << left_contig << "\t=========\t"
            //                   << fix_seq << "\t=========\t"
            //                   << right_contig << std::endl;
            //         {
            //             double indicate_value = CalcMeanAbsoluteContrast(pred_counts, succ_counts);
            //             std::cerr << contig_list.at(pred_serial).GetSequence() << "\t" << contig_list.at(pred_serial).GetTag();
            //             // for (double c : pred_counts)
            //             // {
            //             //     std::cerr << "\t" << c;
            //             // }
            //             std::cerr << std::endl;
            //             std::cerr << contig_list.at(succ_serial).GetSequence() << "\t" << contig_list.at(succ_serial).GetTag();
            //             // for (double c : succ_counts)
            //             // {
            //             //     std::cerr << "\t" << c;
            //             // }
            //             std::cerr << std::endl;
            //             // std::cerr << contig_list.at(pred_serial).GetSequence() << '\t' << contig_list.at(succ_serial).GetSequence();
            //             std::cerr << "\t=========>" << indicate_value << std::endl;
            //         }
            //     }
            // }
            if (interv_method == INTERV_PEARSON && CalcPearsonCorrelation(pred_counts, succ_counts) < interv_thres)
            {
                continue;
            }
            else if (interv_method == INTERV_SPEARMAN && CalcSpearmanCorrelation(pred_counts, succ_counts) < interv_thres)
            {
                continue;
            }
            else if (interv_method == INTERV_MAC && CalcMeanAbsoluteContrast(pred_counts, succ_counts) > interv_thres)
            {
                continue;
            }
        }
        // if (ov_kt.second.IsPredtigRC() && ov_kt.second.IsSucctigRC())
        // {
        //     contig_list.emplace_back(contig_list.at(succ_serial), contig_list.at(pred_serial), n_overlap, quant_mode);
        // }
        // else
        // {
        if (ov_kt.second.IsPredtigRC())
        {
            contig_list.at(pred_serial).ReverseComplement();
        }
        if (ov_kt.second.IsSucctigRC())
        {
            contig_list.at(succ_serial).ReverseComplement();
        }
        contig_list.emplace_back(contig_list.at(pred_serial), contig_list.at(succ_serial), n_overlap, quant_mode);
        // }
        contig_list.at(pred_serial).SetUsed();
        contig_list.at(succ_serial).SetUsed();
        have_new_extensions = true;
    }
    return have_new_extensions;
}

int main(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string kmer_count_path, sample_info_path, interv_method("none"), quant_mode("rep"), tmp_dir("./"), rep_value_cname;
    double interv_thres;
    size_t k_length(31);
    bool stranded(true), disk_mode(false), comp_adj(false);
    size_t min_overlap(15);

    ParseOptions(argc, argv, stranded, disk_mode, k_length, min_overlap, sample_info_path, tmp_dir,
                 interv_method, interv_thres, quant_mode, comp_adj, rep_value_cname, kmer_count_path);
    PrintRunInfo(stranded, disk_mode, k_length, min_overlap, sample_info_path, tmp_dir,
                 interv_method, interv_thres, quant_mode, comp_adj, rep_value_cname, kmer_count_path);

    std::cerr << "Option dealing finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::vector<ContigInfo> contig_list(1, ContigInfo("NOT_A_CONTIG", -1, 0)); // Have a placeholder contig at the index 0
    ColumnInfo column_info;
    std::string count_index_path = tmp_dir + "counts.idx";
    ScanCountTable(contig_list, column_info, kmer_count_path, rep_value_cname, sample_info_path, count_index_path, stranded, disk_mode);

    std::cerr << "Count table Scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::ifstream index_file(count_index_path);
    ExitIf(!index_file.is_open(), "ERROR: k-mer count index file " + count_index_path + " was not found.");
    // PrintContigList(std::cout, contig_list, column_info.GetNbSample(), quant_mode, k_length, stranded, disk_mode, index_file, kmer_disk_pos);

    for (size_t n_overlap(k_length - 1); n_overlap >= min_overlap; --n_overlap)
    {
        std::cerr << "Merging contigs with overlap " << n_overlap << std::endl;
        bool have_new_extensions(true);
        while (have_new_extensions)
        {
            std::cerr << "\tcontig list size: " << contig_list.size() - 1 << std::endl;
            std::map<uint64_t, OverlapKnot> overlap_knot_dict;
            MakeOverlapKnotDict(overlap_knot_dict, contig_list, stranded, n_overlap);
            have_new_extensions = DoExtension(contig_list, overlap_knot_dict, n_overlap, interv_method, interv_thres,
                                              column_info.GetNbSample(), quant_mode, disk_mode, comp_adj, index_file);
            contig_list.erase(std::remove_if(contig_list.begin(), contig_list.end(),
                                             [](ContigInfo contig) { return contig.IsUsed(); }),
                              contig_list.end());
        }
    }
    column_info.PrintSampleNames(std::cout, "contig\tnb_merged_kmers\ttag");
    PrintContigList(std::cout, contig_list, column_info.GetNbSample(), disk_mode, index_file);
    index_file.close();

    std::cerr << "Contig extension finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
