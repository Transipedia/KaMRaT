#ifndef KAMRAT_RUNINFOPARSER_UTILS_HPP
#define KAMRAT_RUNINFOPARSER_UTILS_HPP

#include <set>
#include <string>

const std::set<std::string> EVALU_METHOD_UNIV{"pearson", "spearman", "mac"};
const std::set<std::string> EVALU_MODE_UNIV{"farthest", "worstAdj"};

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
