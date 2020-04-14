#ifndef KAMRAT_QUERY_CONTIGSAMPLEINFO_HPP
#define KAMRAT_QUERY_CONTIGSAMPLEINFO_HPP

#include <string>
#include <unordered_map>
#include <vector>

class ContigSampleInfo
{
public:
    ContigSampleInfo(const std::string &contig_seq);
    const std::string &GetContigSeq() const;
    void AddKmerCounts(const std::string &count_line,
                       const std::vector<std::string> &col_name_vect,
                       const std::vector<bool> &is_col_sample);
    void EstimatePrint(std::ostream &out_s,
                       double (*eval_func)(const std::vector<size_t> &counts),
                       const std::vector<std::string> &col_name_vect,
                       const std::vector<bool> &is_col_sample) const;

private:
    std::string contig_seq_;
    std::unordered_map<std::string, std::vector<size_t>> sample_counts_;
};

#endif //KAMRAT_QUERY_CONTIGSAMPLEINFO_HPP