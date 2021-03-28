#include <iostream>
#include <cstring>

#define RESET "\033[0m"
#define BOLDRED "\033[1m\033[31m"
#define BOLDYELLOW "\033[1m\033[33m"

int IndexMain(int argc, char *argv[]);
// int FilterMain(int argc, char *argv[]);
// int MaskMain(int argc, char *argv[]);
// int MergeMain(int argc, char *argv[]);
// int RankMain(int argc, char *argv[]);
// int EstimateMain(int argc, char *argv[]);

const void Welcome()
{
    std::cerr << "======================================================================" << std::endl
              << " |                                                                  |" << std::endl
              << " |  888    d8P           888b     d888 8888888b.       88888888888  |" << std::endl
              << " |  888   d8P            8888b   d8888 888   Y88b          888      |" << std::endl
              << " |  888  d8P             88888b.d88888 888    888          888      |" << std::endl
              << " |  888d88K      8888b.  888Y88888P888 888   d88P  8888b.  888      |" << std::endl
              << " |  8888888b        \"88b 888 Y888P 888 8888888P\"      \"88b 888      |" << std::endl
              << " |  888  Y88b   .d888888 888  Y8P  888 888 T88b   .d888888 888      |" << std::endl
              << " |  888   Y88b  888  888 888   \"   888 888  T88b  888  888 888      |" << std::endl
              << " |  888    Y88b \"Y888888 888       888 888   T88b \"Y888888 888      |" << std::endl
              << " |                                                                  |" << std::endl
              << "================== k-mer Matrix, Really Tremendous! ==================" << std::endl
              << std::endl;
}

const void PrintHelper()
{
    std::cerr << "[Usage]" << std::endl
              << "    kamrat <command> [options]" << std::endl
              << "[Command]" << std::endl
              << "    index:       index the given matrix, get ready for other modules" << std::endl
              << "    filter:      filter count table according to occurence among samples" << std::endl
              << "    mask:        mask/select count table according to given sequences" << std::endl
              << "    merge:       merge k-mer count table into contig count table" << std::endl
              << "    rank:        rank count table features according to association to sample condition" << std::endl
              << "    estimate:    estimate count for each given sequence" << std::endl
              << std::endl;
}

int main(int argc, char *argv[])
{
    Welcome();
    try
    {
        if (argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "--help") == 0)
        {
            PrintHelper();
        }
        else if (strcmp(argv[1], "index") == 0)
        {
            IndexMain(argc - 1, &(argv[1]));
        }
        // else if (strcmp(argv[1], "filter") == 0)
        // {
        //     FilterMain(argc - 1, &(argv[1]));
        // }
        // else if (strcmp(argv[1], "mask") == 0)
        // {
        //     MaskMain(argc - 1, &(argv[1]));
        // }
        // else if (strcmp(argv[1], "merge") == 0)
        // {
        //     MergeMain(argc - 1, &(argv[1]));
        // }
        // else if (strcmp(argv[1], "rank") == 0)
        // {
        //     RankMain(argc - 1, &(argv[1]));
        // }
        // else if (strcmp(argv[1], "estimate") == 0)
        // {
        //     EstimateMain(argc - 1, &(argv[1]));
        // }
        else
        {
            PrintHelper();
            throw std::invalid_argument("unknown command: " + std::string(argv[1]));
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << BOLDRED << "[ERROR] " << e.what() << RESET << std::endl;
        abort();
    }

    return EXIT_SUCCESS;
}
