#include <iostream>
#include <ctime>
#include <fstream>
#include <map>
#include <unordered_map>
#include <cmath>
#include <algorithm> // std::remove_if

#include "merge/merge_runinfo.hpp"
#include "data_struct/contig_elem.hpp"
#include "data_struct/merge_knot.hpp"

using code2pos_t = std::map<uint64_t, size_t>;
using contigVect_t = std::vector<std::unique_ptr<ContigElem>>;

uint64_t Seq2Int(const std::string &seq, const size_t k_length, const bool stranded); // in utils/seq_coding.cpp
uint64_t GetRC(const uint64_t code, size_t k_length);                                 // in utils/seq_coding.cpp

const double CalcPearsonDist(const std::vector<float> &x, const std::vector<float> &y);  // in utils/vect_opera.cpp
const double CalcSpearmanDist(const std::vector<float> &x, const std::vector<float> &y); // in utils/vect_opera.cpp
const double CalcMACDist(const std::vector<float> &x, const std::vector<float> &y);      // in utils/vect_opera.cpp

void LoadIndexMeta(size_t &nb_smp, size_t &k_len, bool &stranded,
                   std::vector<std::string> &colname_vect, std::vector<double> &smp_sum_vect,
                   const std::string &idx_meta_path);                             // in utils/index_loading.cpp
void LoadPosVect(std::vector<size_t> &pos_vect, const std::string &idx_pos_path); // in utils/index_loading.cpp
const std::vector<float> &GetCountVect(std::vector<float> &count_vect,
                                       std::ifstream &idx_mat, const size_t pos, const size_t nb_smp); // in utils/index_loading.cpp

// ----- when the -select argument is provided: merging on subset of index ----- //
void LoadSelect(std::unordered_map<uint64_t, float> &sel_kmer_map, const std::string &sel_path, const std::string &rep_mode,
                const size_t k_len, const bool stranded)
{
    std::ifstream sel_file(sel_path);
    if (!sel_file.is_open())
    {
        throw std::invalid_argument("path for input sequences not accessable: " + sel_path);
    }
    std::string line;
    size_t split_pos;
    float rep_val(0);
    while (std::getline(sel_file, line))
    {
        split_pos = line.find_first_of("\t ");
        if (split_pos != std::string::npos)
        {
            rep_val = std::stof(line.substr(split_pos + 1));
            line = line.substr(0, split_pos);
            if (rep_mode == "minabs")
            {
                rep_val = fabs(rep_val);
            }
            else if (rep_mode == "max")
            {
                rep_val = -rep_val;
            }
            else if (rep_mode == "maxabs")
            {
                rep_val = -fabs(rep_val);
            }
        }
        if (!sel_kmer_map.insert({Seq2Int(line, k_len, stranded), rep_val}).second) // check unicity of given k-mers
        {
            throw std::domain_error("checking unicity failed: k-mer not unique in the given sequence list: " + line);
        }
    }
    sel_file.close();
}

void LoadKMerPos(std::map<uint64_t, size_t> &)

void InitializeContigList(contigVect_t &ctg_vect, std::ifstream &idx_mat, const std::string &idx_pos_path,
                          const std::string &sel_path, const std::string &rep_mode, const size_t k_len, const bool stranded)
{
    std::unordered_map<uint64_t, float> sel_kmer_map;
    LoadSelect(sel_kmer_map, sel_path, rep_mode, k_len, stranded);
    size_t pos;
    while (idx_pos.read(reinterpret_cast<char *>(&pos), sizeof(size_t)))
    {
        pos_vect.emplace_back(pos);
    }
    std::string term;
    for (size_t p : pos_vect) // erase while loading
    {
        idx_mat.ignore(nb_smp * sizeof(float)); // skip the indexed count vector
        idx_mat >> term;
        ctg_vect.emplace_back(std::make_unique<ContigElem>(term, p, 0));
    }
    idx_pos.close();
}

// ----- load all k-mers in index in the order of uint64 code ----- //
void MakeKMerMap(code2pos_t &code_pos_map, std::ifstream &idx_pos, std::ifstream &idx_mat,
                 const size_t k_len, const bool stranded, const size_t nb_smp)
{
    std::string kmer;
    size_t pos;
    while (idx_pos.read(reinterpret_cast<char *>(&pos), sizeof(size_t)))
    {
        idx_mat.seekg(pos + sizeof(float) * nb_smp);
        idx_mat >> kmer;
        if (!code_pos_map.insert({Seq2Int(kmer, k_len, stranded), pos}).second)
        {
            throw std::domain_error("checking unicity failed: " + kmer + " is not unique in the given sequence list");
        }
    }
    std::vector<float> count_vect;
    for (const auto &elem : code_pos_map)
    {
        std::cout << elem.first;
        for (const float x : GetCountVect(count_vect, idx_mat, elem.second, nb_smp))
        {
            std::cout << "\t" << x;
        }
        std::cout << std::endl;
    }
}


// ----- when the -inseq argument is not provided: merging on all index ----- //
void MakeContigFromIndex(contigVect_t &ctg_vect, std::ifstream &idx_mat, const std::vector<size_t> &pos_vect, size_t nb_smp)
{
    std::string term;
    for (size_t p : pos_vect) // erase while loading
    {
        idx_mat.ignore(nb_smp * sizeof(float)); // skip the indexed count vector
        idx_mat >> term;
        ctg_vect.emplace_back(std::make_unique<ContigElem>(term, p, 0));
    }
}

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
        hashed_merge_knots.insert({prefix, MergeKnot()}).first->second.AddContig(i_ctg, is_prefix_rc, (is_prefix_rc ? "pred" : "succ"));
        hashed_merge_knots.insert({suffix, MergeKnot()}).first->second.AddContig(i_ctg, is_suffix_rc, (is_suffix_rc ? "succ" : "pred"));
    }
}

const bool DoExtension(contigVect_t &ctg_vect, fix2knot_t &hashed_mergeknot_list, const size_t i_ovlp,
                       const std::string &interv_method, const float interv_thres,
                       std::ifstream &idx_mat, const size_t nb_smp)
{
    static std::vector<float> pred_counts, succ_counts;
    bool has_new_extensions(false);
    for (auto it = hashed_mergeknot_list.cbegin(); it != hashed_mergeknot_list.cend(); ++it)
    {
        if (!it->second.IsMergeable())
        {
            continue;
        }
        if (ctg_vect[it->second.GetSerial("pred")] == nullptr || ctg_vect[it->second.GetSerial("succ")] == nullptr)
        {
            continue;
        }
        bool pred_rc = it->second.IsRC("pred"), succ_rc = it->second.IsRC("succ");
        if (interv_method == "pearson" &&
            CalcPearsonDist(GetCountVect(pred_counts, idx_mat, ctg_vect[it->second.GetSerial("pred")]->GetRearPos(pred_rc), nb_smp),
                            GetCountVect(succ_counts, idx_mat, ctg_vect[it->second.GetSerial("succ")]->GetHeadPos(succ_rc), nb_smp)) >= interv_thres)
        {
            continue;
        }
        else if (interv_method == "spearman" &&
                 CalcSpearmanDist(GetCountVect(pred_counts, idx_mat, ctg_vect[it->second.GetSerial("pred")]->GetRearPos(pred_rc), nb_smp),
                                  GetCountVect(succ_counts, idx_mat, ctg_vect[it->second.GetSerial("succ")]->GetHeadPos(succ_rc), nb_smp)) >= interv_thres)
        {
            continue;
        }
        else if (interv_method == "mac" &&
                 CalcMACDist(GetCountVect(pred_counts, idx_mat, ctg_vect[it->second.GetSerial("pred")]->GetRearPos(pred_rc), nb_smp),
                             GetCountVect(succ_counts, idx_mat, ctg_vect[it->second.GetSerial("succ")]->GetHeadPos(succ_rc), nb_smp)) >= interv_thres)
        {
            continue;
        }
        if (pred_rc && succ_rc)
        {
            ctg_vect[it->second.GetSerial("succ")]->RightExtend(std::move(ctg_vect[it->second.GetSerial("pred")]), false, i_ovlp);
        }
        else if (pred_rc)
        {
            ctg_vect[it->second.GetSerial("succ")]->LeftExtend(std::move(ctg_vect[it->second.GetSerial("pred")]), true, i_ovlp);
        }
        else
        {
            ctg_vect[it->second.GetSerial("pred")]->RightExtend(std::move(ctg_vect[it->second.GetSerial("succ")]), succ_rc, i_ovlp);
        }
        // it = hashed_mergeknot_list.erase(it);
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

void PrintContigVect(const contigVect_t &ctg_vect, std::ifstream &idx_mat, const size_t nb_smp)
{
    std::vector<float> count_vect;
    for (const auto &elem : ctg_vect)
    {
        std::cout << elem->GetSeq();
        for (const float x : GetCountVect(count_vect, idx_mat, elem->GetHeadPos(false), nb_smp))
        {
            std::cout << "\t" << x;
        }
        std::cout << std::endl;
        // std::cout << elem->GetSeq();
        // for (const float x : GetCountVect(count_vect, idx_mat, elem->GetRearPos(false), nb_smp))
        // {
        //     std::cout << "\t" << x;
        // }
        // std::cout << std::endl;
    }
}

int MergeMain(int argc, char **argv)
{
    MergeWelcome();

    std::clock_t begin_time = clock(), inter_time;
    std::string idx_dir, sel_path, rep_mode("min"), design_path, itv_mthd("none"), out_path;
    float itv_thres(0.5);
    size_t max_ovlp(0), min_ovlp(0);

    ParseOptions(argc, argv, idx_dir, max_ovlp, min_ovlp, sel_path, rep_mode, itv_mthd, itv_thres, out_path);

    size_t nb_smp, k_len;
    bool stranded;
    std::vector<std::string> colname_vect;
    {
        std::vector<double> smp_sum_vect; // smp_sum_vect not needed in KaMRaT-merge
        LoadIndexMeta(nb_smp, k_len, stranded, colname_vect, smp_sum_vect, idx_dir + "/idx-meta.bin");
        if (k_len == 0)
        {
            throw std::domain_error("KaMRaT-merge relies on the index in k-mer mode, please rerun KaMRaT-index with -klen option");
        }
        if (k_len <= max_ovlp)
        {
            throw std::invalid_argument("max overlap (" + std::to_string(max_ovlp) + ") should not exceed k-mer length (" + std::to_string(k_len) + ")");
        }
        PrintRunInfo(idx_dir, k_len, max_ovlp, min_ovlp, stranded, sel_path, rep_mode, itv_mthd, itv_thres, out_path);
    }
    std::vector<size_t> pos_vect;
    LoadPosVect(pos_vect, idx_dir + "/idx-pos.bin");
    // std::cout << nb_smp << std::endl;
    // std::vector<bool> is_smp_col(nb_smp, true);
    // if (!design_path.empty())
    // {
    //     LoadDesign(is_smp_col, design_path, colname_vect);
    // }
    // for (bool x : is_smp_col)
    // {
    //     std::cout << (x ? "T" : "F") << "\t";
    // }
    // std::cout << std::endl;

    // code2kmer_t code_kmer_map;
    // MakeKMerMap(code_kmer_map, idx_dir + "/idx-pos.bin", k_len, stranded);

    std::cerr << "Option parsing and index loading finished, execution time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::ifstream idx_pos(idx_dir + "/idx-pos.bin"), idx_mat(idx_dir + "/idx-mat.bin");
    if (!idx_pos.is_open())
    {
        throw std::invalid_argument("loading idx-pos failed, KaMRaT index folder not found or may be corrupted");
    }
    if (!idx_mat.is_open())
    {
        throw std::invalid_argument("loading index-mat failed, KaMRaT index folder not found or may be corrupted");
    }

    code2pos_t code_pos_map;
    MakeKMerMap(code_pos_map, idx_pos, idx_mat, k_len, stranded, nb_smp);
    if (!sel_path.empty())
    {
        // LoadSelectList();
    }

    contigVect_t ctg_vect;
    if (sel_path.empty())
    {
        MakeContigFromIndex(ctg_vect, idx_mat, pos_vect, nb_smp);
    }
    else
    {
        // MakeContigFromList(ctg_vect, sel_path, code_kmer_map, k_len, stranded);
        // code_kmer_map.clear();
    }

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
            has_new_extensions = DoExtension(ctg_vect, hashed_merge_knots, i_ovlp, itv_mthd, itv_thres, idx_mat, nb_smp);
            // PrintContigList(ctg_vect, kmer_count_tab, k_len, stranded, min_ovl, interv_method, interv_thres, quant_mode, idx_mat);
            ctg_vect.erase(std::remove_if(ctg_vect.begin(), ctg_vect.end(), [](const auto &elem) { return elem == nullptr; }),
                           ctg_vect.end());
            hashed_merge_knots.clear();
        }
    }
    PrintContigVect(ctg_vect, idx_mat, nb_smp);
    idx_mat.close();

    std::cerr << "Contig print finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}