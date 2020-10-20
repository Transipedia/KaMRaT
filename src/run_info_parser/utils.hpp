#ifndef KAMRAT_RUNINFOPARSER_UTILS_HPP
#define KAMRAT_RUNINFOPARSER_UTILS_HPP

#include <set>
#include <string>

#define MAX_PEARSON_DEFAULT 0.20
#define MAX_SPEARMAN_DEFAULT 0.26
#define MAX_MAC_DEFAULT 0.25

const std::set<std::string> EVALU_METHOD_UNIV{"pearson", "spearman", "mac"};
const std::set<std::string> EVALU_MODE_UNIV{"farthest", "worstAdj"};
const std::set<std::string> INTERV_METHOD_UNIV{"pearson", "spearman", "mac"};
const std::set<std::string> QUANT_MODE_UNIV{"rep", "mean"};
const std::set<std::string> SORT_MODE_UNIV{"dec", "inc", "dec:abs", "inc:abs"};

inline void SubCommandParser(std::string &command, std::string &sub_command)
{
    size_t split_pos = command.find(":");
    if (split_pos != std::string::npos)
    {
        sub_command = command.substr(split_pos + 1);
        command = command.substr(0, split_pos);
    }
}

#endif //KAMRAT_RUNINFOPARSER_UTILS_HPP
