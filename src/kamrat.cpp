#include <iostream>
#include <cstring>

int FilterMain(int argc, char *argv[]);
int MaskMain(int argc, char *argv[]);
int MergeMain(int argc, char *argv[]);
int RankMain(int argc, char *argv[]);

const void Welcome()
{
    std::cerr << "===================================================================================" << std::endl
              << "******                           Welcome to KaMRaT                           ******" << std::endl
              << "******                                                                       ******" << std::endl
              << "****** --------------------------------------------------------------------- ******" << std::endl
	      << std::endl;
}

const void PrintHelper()
{
    std::cerr << "[Usage]" << std::endl 
              << "    kamrat <command> [options] COUNT_TAB_PATH" << std::endl
              << "[Command]" << std::endl
              << "    filter:    filter count table according to occurence among samples" << std::endl
              << "    mask:      mask/select count table according to given sequences" << std::endl
              << "    merge:     merge k-mer count table into contig count table" << std::endl
              << "    rank:      rank count table features according to association to sample condition" << std::endl
              << std::endl;
}

int main(int argc, char *argv[])
{
    Welcome();
    if (argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "--help") == 0)
    {
        PrintHelper();
    }
    else if (strcmp(argv[1], "filter") == 0) {
        FilterMain(argc - 1, &(argv[1]));
    }
    else if (strcmp(argv[1], "mask") == 0) {
        MaskMain(argc - 1, &(argv[1]));
    }
    else if (strcmp(argv[1], "merge") == 0) {
        MergeMain(argc - 1, &(argv[1]));
    }
    else if (strcmp(argv[1], "rank") == 0) {
        RankMain(argc - 1, &(argv[1]));
    }
    else
    {
        PrintHelper();
        throw std::invalid_argument("unknown command: " + std::string(argv[1]));
    }
    return EXIT_SUCCESS;
}
