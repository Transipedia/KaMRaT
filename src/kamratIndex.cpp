#include <iostream>
#include <ctime>
#include <string>

#include "index/index_runinfo.hpp"

void ScanIndex(const std::string &count_tab_path)
{
}

int MergeMain(int argc, char **argv)
{
    std::clock_t begin_time = clock(), inter_time;
    std::string out_dir, count_tab_path;

    ParseOptions(argc, argv, out_dir, count_tab_path);
    PrintRunInfo(out_dir, count_tab_path);

    std::cerr << "Count table scanning finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    inter_time = clock();

    std::cerr << "Contig print finished, execution time: " << (float)(clock() - inter_time) / CLOCKS_PER_SEC << "s." << std::endl;
    std::cerr << "Total executing time: " << (float)(clock() - begin_time) / CLOCKS_PER_SEC << "s." << std::endl;

    return EXIT_SUCCESS;
}
