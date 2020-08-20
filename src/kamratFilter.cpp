#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <ctime>

#include "data_struct/count_tab_header.hpp"
#include "run_info_parser/filter.hpp"

const long GetSampleSerial(const std::string &sample_name, const CountTabHeader &count_tab_header)
{
    for (size_t i(0); i < count_tab_header.GetNbColumn(); ++i)
    {
        if (count_tab_header.GetColName(i) == sample_name)
        {
            return count_tab_header.GetColSerial(i);
        }
    }
    throw std::domain_error("unknown sample names for expressed sample filter: " + sample_name);
}

const bool IsSampleToCount(const std::string &this_level,
                           const std::string &this_name,
                           const std::string &other_level,
                           const std::string &other_name,
                           const CountTabHeader &count_tab_header,
                           const size_t i_col)
{
    if (this_name == "all") // all samples
    {
        return true;
    }
    else if (this_name == "rest" && other_level == "smp" && !count_tab_header.CheckColName(i_col, other_name)) // samples other than silent sample
    {
        return true;
    }
    else if (this_name == "rest" && other_level == "cond" && !count_tab_header.CheckColCondition(i_col, other_name)) // simples other than silent condition
    {
        return true;
    }
    else if (this_level == "smp" && count_tab_header.CheckColName(i_col, this_name)) // expressed sample
    {
        return true;
    }
    else if (this_level == "cond" && count_tab_header.CheckColCondition(i_col, this_name)) // samples in expressed condition
    {
        return true;
    }
    return false;
}

void ScanCountTab(const std::string &count_tab_path,
                  const std::string &sample_info_path,
                  const std::string &express_level,
                  const std::string &express_name,
                  const size_t express_min_abd,
                  const size_t express_min_rec,
                  const std::string &silent_level,
                  const std::string &silent_name,
                  const size_t silent_max_abd,
                  const size_t silent_min_rec)
{
    std::ifstream count_tab_file(count_tab_path);
    if (!count_tab_file.is_open())
    {
        throw std::domain_error("count table file " + count_tab_path + " is not found");
    }
    CountTabHeader count_tab_header;
    //----- Dealing with Header Line for Constructing ColumnInfo Object -----//
    std::set<std::string> preserved_condition_names{"all", "rest"};
    std::string line;
    std::getline(count_tab_file, line);
    count_tab_header.MakeSmpCond(sample_info_path, preserved_condition_names);
    count_tab_header.MakeColumnInfo(line, "??%#NOT_CARE_SCORE-@_@");
    const long express_smp_serial = (express_level == "smp" ? GetSampleSerial(express_name, count_tab_header) : -1),
               silent_smp_serial = (silent_level == "smp" ? GetSampleSerial(silent_level, count_tab_header) : -1);
    std::vector<size_t> label_vect;
    count_tab_header.GetSmpLabels(label_vect);
    std::cout << line << std::endl;
    //----- Dealing with Following k-mer Count Lines -----//
    while (std::getline(count_tab_file, line))
    {
        std::istringstream conv(line);
        std::string term;
        size_t express_rec(0), silent_rec(0);
        for (size_t i(0); conv >> term && i < count_tab_header.GetNbColumn(); ++i)
        {
            if (!count_tab_header.IsSample(i))
            {
                continue;
            }
            if (IsSampleToCount(express_level, express_name, silent_level, silent_name, count_tab_header, i) && std::stof(term) >= express_min_abd)
            {
                ++express_rec;
            }
            if (IsSampleToCount(silent_level, silent_name, express_level, express_name, count_tab_header, i) && std::stof(term) <= silent_max_abd)
            {
                ++silent_rec;
            }
            if (conv.fail() || i == count_tab_header.GetNbColumn())
            {
                throw std::domain_error("line parsing not coherent with header");
            }
        }
        // std::cerr << express_rec << "\t" << silent_rec << std::endl;
        if (express_rec >= express_min_rec && silent_rec >= silent_min_rec)
        {
            std::cout << line << std::endl;
        }
    }
    count_tab_file.close();
}

int main(int argc, char *argv[])
{
    std::clock_t begin_time = clock();
    std::string count_tab_path, sample_info_path, express_level, express_name, silent_level, silent_name;
    int express_min_abd(-1), express_min_rec(-1), silent_max_abd(-1), silent_min_rec(-1);

    ParseOptions(argc, argv, count_tab_path, express_level, express_name, express_min_rec, express_min_abd, silent_level, silent_name, silent_min_rec, silent_max_abd, sample_info_path);
    PrintRunInfo(count_tab_path, express_level, express_name, express_min_rec, express_min_abd, silent_level, silent_name, silent_min_rec, silent_max_abd, sample_info_path);

    ScanCountTab(count_tab_path, sample_info_path, express_level, express_name, express_min_abd, express_min_rec, silent_level, silent_name, silent_max_abd, silent_min_rec);

    std::cerr << "Executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
    return EXIT_SUCCESS;
}
