#include <iostream>
#include <ctime>
#include <fstream>
#include <map>
#include <unordered_map>
#include <cmath>
#include <algorithm> // std::remove_if

#include "merge_runinfo.hpp"
#include "contig_elem.hpp"
#include "merge_knot.hpp"
#include "seq_coding.hpp"
#include "index_loading.hpp"

using contigVect_t = std::vector<std::unique_ptr<ContigElem>>;

const double CalcPearsonDist(const std::vector<float> &x, const std::vector<float> &y);  // in utils/vect_opera.cpp
const double CalcSpearmanDist(const std::vector<float> &x, const std::vector<float> &y); // in utils/vect_opera.cpp
const double CalcMACDist(const std::vector<float> &x, const std::vector<float> &y);      // in utils/vect_opera.cpp

const bool MakeContigListFromIndex(contigVect_t &ctg_vect, const std::string &idx_pos_path,
                                   std::ifstream &idx_mat, const size_t nb_smp)
{
    std::ifstream idx_pos(idx_pos_path);
    if (!idx_pos.is_open())
    {
        throw std::invalid_argument("loading index-pos failed, KaMRaT index folder not found or may be corrupted");
    }
    size_t rep_pos;
    std::string kmer_seq;
    while (idx_pos.read(reinterpret_cast<char *>(&rep_pos), sizeof(uint64_t)) && idx_pos.read(reinterpret_cast<char *>(&rep_pos), sizeof(size_t)))
    {
        idx_mat.seekg(rep_pos + nb_smp * sizeof(float)); // skip the indexed count vector
        idx_mat >> kmer_seq;
        ctg_vect.emplace_back(std::make_unique<ContigElem>(kmer_seq, rep_pos, 0));
    }
    idx_pos.close();
    return false;
}

const bool MakeContigListFromFile(contigVect_t &ctg_vect, const std::string &with_path)
{
    std::ifstream with_file(with_path);
    if (!with_file.is_open())
    {
        throw std::invalid_argument("path for input sequences not accessable: " + with_path);
    }
    std::string kmer_seq;
    size_t nb_kmer, rep_pos;
    float rep_val(0);
    bool has_value(false), read_succ(false);
    while (with_file >> kmer_seq >> rep_val >> nb_kmer)
    {
        read_succ = true;
        with_file.ignore(1);
        with_file.read(reinterpret_cast<char *>(&rep_pos), sizeof(size_t));
        if (nb_kmer > 1)
        {
            throw std::domain_error("find a feature with k-mer number > 1: " + kmer_seq);
        }
        if (rep_val > 1E-5 || rep_val < -1E-5)
        {
            has_value = true;
        }
        ctg_vect.emplace_back(std::make_unique<ContigElem>(kmer_seq, rep_pos, rep_val));
    }
    with_file.close();
    if (!read_succ)
    {
        throw std::invalid_argument("not valid file for input sequences: " + with_path);
    }
    return has_value;
}

/**
 */
void MakeOverlapKnots(fix2knot_t &hashed_merge_knots, const contigVect_t &ctg_vect, const bool stranded, const size_t i_ovlp)
{
    for (size_t i_ctg(0); i_ctg < ctg_vect.size(); ++i_ctg)
    {
        const std::string &seq = ctg_vect[i_ctg]->GetSeq();
        uint64_t prefix = Seq2Int(seq, i_ovlp, true), suffix = Seq2Int(seq.substr(seq.size() - i_ovlp), i_ovlp, true);
        bool is_prefix_rc(false), is_suffix_rc(false);
        if (!stranded)
        {
            uint64_t prefix_rc = GetRC(prefix, i_ovlp), suffix_rc = GetRC(suffix, i_ovlp);
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
        if (prefix == suffix) // if the k-mer has equal prefix and suffix
        {
            const auto &ins_pair = hashed_merge_knots.insert({prefix, MergeKnot()});
            if (!ins_pair.first->second.HasPred() && !is_suffix_rc)
            {
                ins_pair.first->second.AddContig(i_ctg, false, "pred");
            }
            else if (!ins_pair.first->second.HasPred() && is_prefix_rc)
            {
                ins_pair.first->second.AddContig(i_ctg, true, "pred");
            }
            else if (!ins_pair.first->second.HasSucc() && !is_prefix_rc)
            {
                ins_pair.first->second.AddContig(i_ctg, false, "succ");
            }
            else if (!ins_pair.first->second.HasSucc() && is_suffix_rc)
            {
                ins_pair.first->second.AddContig(i_ctg, true, "succ");
            }
            else
            {
                ins_pair.first->second.AddContig(i_ctg, is_prefix_rc, (is_prefix_rc ? "pred" : "succ"));
            }
        }
        else // if the k-mer has different prefix and suffix
        {
            hashed_merge_knots.insert({prefix, MergeKnot()}).first->second.AddContig(i_ctg, is_prefix_rc, (is_prefix_rc ? "pred" : "succ"));
            hashed_merge_knots.insert({suffix, MergeKnot()}).first->second.AddContig(i_ctg, is_suffix_rc, (is_suffix_rc ? "succ" : "pred"));
        }
    }
}

const bool IsFirstContigRep(const float ctg_val1, const float ctg_val2, const std::string &rep_mode)
{
    if (rep_mode == "min")
    {
        return (ctg_val1 <= ctg_val2);
    }
    else if (rep_mode == "max")
    {
        return (ctg_val1 >= ctg_val2);
    }
    else if (rep_mode == "minabs")
    {
        return (fabs(ctg_val1) <= fabs(ctg_val2));
    }
    else // maxabs
    {
        return (fabs(ctg_val1) >= fabs(ctg_val2));
    }
}

const bool DoExtension(contigVect_t &ctg_vect, const fix2knot_t &hashed_mergeknot_list, const size_t i_ovlp,
                       const std::string &interv_method, const float interv_thres,
                       std::ifstream &idx_mat, const size_t nb_smp, const std::string &rep_mode)
{
    static std::vector<float> pred_counts, succ_counts;
    bool has_new_extensions(false);
    for (auto it = hashed_mergeknot_list.cbegin(); it != hashed_mergeknot_list.cend(); ++it)
    {
        if (!it->second.IsMergeable())
        {
            continue;
        }
        auto &pred_ctg = ctg_vect[it->second.GetSerial("pred")], &succ_ctg = ctg_vect[it->second.GetSerial("succ")];
        if (pred_ctg == nullptr || succ_ctg == nullptr)
        {
            continue;
        }
        const bool pred_rc = it->second.IsRC("pred"), succ_rc = it->second.IsRC("succ");
        if (interv_method == "pearson" &&
            CalcPearsonDist(GetCountVect(pred_counts, idx_mat, pred_ctg->GetRearPos(pred_rc), nb_smp),
                            GetCountVect(succ_counts, idx_mat, succ_ctg->GetHeadPos(succ_rc), nb_smp)) >= interv_thres)
        {
            continue;
        }
        else if (interv_method == "spearman" &&
                 CalcSpearmanDist(GetCountVect(pred_counts, idx_mat, pred_ctg->GetRearPos(pred_rc), nb_smp),
                                  GetCountVect(succ_counts, idx_mat, succ_ctg->GetHeadPos(succ_rc), nb_smp)) >= interv_thres)
        {
            continue;
        }
        else if (interv_method == "mac" &&
                 CalcMACDist(GetCountVect(pred_counts, idx_mat, pred_ctg->GetRearPos(pred_rc), nb_smp),
                             GetCountVect(succ_counts, idx_mat, succ_ctg->GetHeadPos(succ_rc), nb_smp)) >= interv_thres)
        {
            continue;
        }
        // the base contig should have minimum p-value or input order //
        if (IsFirstContigRep(pred_ctg->GetRepVal(), succ_ctg->GetRepVal(), rep_mode)) // merge right to left
        {
            if (pred_rc) // prevent base contig from reverse-complement transformation, for being coherent with merging knot
            {
                pred_ctg->LeftExtend(std::move(succ_ctg), !succ_rc, i_ovlp); // reverse this ctg and extend at right <=> reverse right ctg and extend to left of this ctg
            }
            else
            {
                pred_ctg->RightExtend(std::move(succ_ctg), succ_rc, i_ovlp);
            }
        }
        else // merge left to right
        {
            if (succ_rc) // prevent base contig from reverse-complement transformation, for being coherent with merging knot
            {
                succ_ctg->RightExtend(std::move(pred_ctg), !pred_rc, i_ovlp); // reverse this ctg and extend at left <=> reverse left ctg and extend to right of this ctg
            }
            else
            {
                succ_ctg->LeftExtend(std::move(pred_ctg), pred_rc, i_ovlp);
            }
        }
        has_new_extensions = true;
    }
    return has_new_extensions;
}

void Int2Seq(std::string &seq, const uint64_t code, const size_t k_length);
void PrintMergeKnots(const fix2knot_t &hashed_merge_knots, const contigVect_t &ctg_vect, const size_t k_len)
{
    for (const auto &elem : hashed_merge_knots)
    {
        if (elem.second.IsMergeable())
        {
            std::string fix,
                contig_pred = ctg_vect[elem.second.GetSerial("pred")]->GetSeq(),
                contig_succ = ctg_vect[elem.second.GetSerial("succ")]->GetSeq();
            Int2Seq(fix, elem.first, k_len);
            std::cout << fix << ": " << contig_pred << " ======= " << contig_succ << std::endl;
        }
    }
}

void PrintHeader(const bool has_value, const std::vector<std::string> &colname_vect)
{
    std::cout << "contig\tnb-merged-kmer";
    if (has_value)
    {
        std::cout << "\trep-value";
    }
    for (const auto &s : colname_vect)
    {
        std::cout << "\t" << s;
    }
    std::cout << std::endl;
}

void PrintWithCounts(const bool has_value, const contigVect_t &ctg_vect, std::ifstream &idx_mat,
                     const std::string &out_mode, const size_t nb_smp, const size_t min_nbkmer)
{
    std::string rep_seq;
    std::vector<float> count_vect;
    for (const auto &elem : ctg_vect)
    {
        if (elem->GetNbMemKmer() < min_nbkmer)
        {
            continue;
        }
        std::cout << elem->GetSeq();
        std::cout << "\t" << elem->GetNbMemKmer();
        if (has_value)
        {
            std::cout << "\t" << elem->GetRepVal();
        }
        std::cout << "\t" << GetTagSeq(rep_seq, idx_mat, elem->GetRepPos(), nb_smp);
        if (out_mode == "rep")
        {
            GetCountVect(count_vect, idx_mat, elem->GetRepPos(), nb_smp); // output sample count vector of representative k-mer
        }
        else if (out_mode == "mean")
        {
            GetMeanCountVect(count_vect, idx_mat, nb_smp, elem->GetMemPosVect()); // output mean sample count vector
        }
        else // out_mode == "median"
        {
            GetMedianCountVect(count_vect, idx_mat, nb_smp, elem->GetMemPosVect()); // output median sample count vector
        }
        for (const float x : count_vect)
        {
            std::cout << "\t" << x;
        }
        std::cout << std::endl;
        count_vect.clear();
    }
}

void PrintAsIntermediate(const contigVect_t &ctg_vect, const size_t min_nbkmer)
{
    size_t rep_pos;
    for (const auto &elem : ctg_vect)
    {
        if (elem->GetNbMemKmer() < min_nbkmer)
        {
            continue;
        }
        std::cout << elem->GetSeq() << "\t0\t" << elem->GetNbMemKmer() << "\t";
        rep_pos = elem->GetRepPos();
        std::cout.write(reinterpret_cast<char *>(&rep_pos), sizeof(size_t));
        for (size_t p : elem->GetMemPosVect())
        {
            if (p != rep_pos) // not repeat the representative position
            {
                std::cout.write(reinterpret_cast<char *>(&p), sizeof(size_t));
            }
        }
        std::cout << std::endl;
    }
}

int MergeMain(int argc, char **argv)
{
    MergeWelcome();

    std::clock_t begin_time = clock(), inter_time;
    std::string idx_dir, with_path, rep_mode("min"), itv_mthd("pearson"), out_path, out_mode;
    float itv_thres(0.20);
    size_t max_ovlp(0), min_ovlp(0), nb_smp(0), k_len(0), min_nbkmer(1);
    bool stranded(false);
    std::vector<std::string> colname_vect;
    ParseOptions(argc, argv, idx_dir, max_ovlp, min_ovlp, with_path, rep_mode, itv_mthd, itv_thres, min_nbkmer, out_path, out_mode);

    // --- Loading ---
    LoadIndexMeta(nb_smp, k_len, stranded, colname_vect, idx_dir + "/idx-meta.bin");
    if (k_len == 0)
    {
        throw std::invalid_argument("KaMRaT-merge relies on the index in k-mer mode, please rerun KaMRaT-index with -klen option");
    }
    if (max_ovlp == 0 && min_ovlp == 0) 
    {
        max_ovlp = k_len - 1;
        min_ovlp = static_cast<size_t>(k_len / 2);
    }
    if (k_len <= max_ovlp)
    {
        throw std::invalid_argument("max overlap (" + std::to_string(max_ovlp) + ") should not exceed k-mer length (" + std::to_string(k_len) + ")");
    }
    PrintRunInfo(idx_dir, k_len, stranded, max_ovlp, min_ovlp, with_path, rep_mode, itv_mthd, itv_thres, min_nbkmer, out_path, out_mode);

    contigVect_t ctg_vect;
    std::ifstream idx_mat(idx_dir + "/idx-mat.bin");
    if (!idx_mat.is_open())
    {
        throw std::invalid_argument("loading index-mat failed, KaMRaT index folder not found or may be corrupted");
    }
    const bool has_value = (with_path.empty() ? MakeContigListFromIndex(ctg_vect, idx_dir + "/idx-pos.bin", idx_mat, nb_smp)
                                              : MakeContigListFromFile(ctg_vect, with_path));
    std::cerr << "Option parsing and index loading finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    fix2knot_t hashed_merge_knots;
    for (size_t i_ovlp(max_ovlp); i_ovlp >= min_ovlp; --i_ovlp)
    {
        std::cerr << "Merging contigs with overlap " << i_ovlp << std::endl;
        bool has_new_extensions(true);
        while (has_new_extensions)
        {
            std::cerr << "\tcontig list size: " << ctg_vect.size() << std::endl;
            MakeOverlapKnots(hashed_merge_knots, ctg_vect, stranded, i_ovlp);
            // PrintMergeKnots(hashed_merge_knots, ctg_vect, k_len);
            has_new_extensions = DoExtension(ctg_vect, hashed_merge_knots, i_ovlp, itv_mthd, itv_thres, idx_mat, nb_smp, rep_mode);
            // PrintContigList(ctg_vect, kmer_count_tab, k_len, stranded, min_ovl, interv_method, interv_thres, quant_mode, idx_mat);
            ctg_vect.erase(std::remove_if(ctg_vect.begin(), ctg_vect.end(), [](const auto &elem)
                                          { return elem == nullptr; }),
                           ctg_vect.end());
            hashed_merge_knots.clear();
        }
    }
    std::cerr << "Contig extension finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

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
    if (has_value && out_mode.empty())
    {
        throw std::invalid_argument("output as intermediate after rank-merge is not possible");
    }
    if (!out_mode.empty())
    {
        PrintHeader(has_value, colname_vect);
        PrintWithCounts(has_value, ctg_vect, idx_mat, out_mode, nb_smp, min_nbkmer);
    }
    else
    {
        PrintAsIntermediate(ctg_vect, min_nbkmer);
    }
    idx_mat.close();

    std::cout.rdbuf(backup_buf);
    if (out_file.is_open())
    {
        out_file.close();
    }
    std::cerr << "Contig print finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    return EXIT_SUCCESS;
}