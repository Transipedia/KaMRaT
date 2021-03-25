#include <iostream>
#include <cstring>

int IndexMain(int argc, char *argv[]);
int FilterMain(int argc, char *argv[]);
int MaskMain(int argc, char *argv[]);
int MergeMain(int argc, char *argv[]);
int RankMain(int argc, char *argv[]);
int EstimateMain(int argc, char *argv[]);

const void Welcome()
{
    std::cerr << "     __    __            __       __  _______           ________     " << std::endl
              << "    /  |  /  |          /  \     /  |/       \         /        |    " << std::endl
              << "    $$ | /$$/   ______  $$  \   /$$ |$$$$$$$  |  ______$$$$$$$$/     " << std::endl
              << "    $$ |/$$/   /      \ $$$  \ /$$$ |$$ |__$$ | /      \  $$ |       " << std::endl
              << "    $$  $$<    $$$$$$  |$$$$  /$$$$ |$$    $$<  $$$$$$  | $$ |       " << std::endl
              << "    $$$$$  \   /    $$ |$$ $$ $$/$$ |$$$$$$$  | /    $$ | $$ |       " << std::endl
              << "    $$ |$$  \ /$$$$$$$ |$$ |$$$/ $$ |$$ |  $$ |/$$$$$$$ | $$ |       " << std::endl
              << "    $$ | $$  |$$    $$ |$$ | $/  $$ |$$ |  $$ |$$    $$ | $$ |       " << std::endl
              << "    $$/   $$/  $$$$$$$/ $$/      $$/ $$/   $$/  $$$$$$$/  $$/        " << std::endl
              << "                                                                     " << std::endl
              << " ================ k-mer matrix, really tremendous ! ================ " << std::endl;
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
    if (argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "--help") == 0)
    {
        PrintHelper();
    }
    else if (strcmp(argv[1], "index") == 0)
    {
        IndexMain(argc - 1, &(argv[2]));
    }
    else if (strcmp(argv[1], "filter") == 0)
    {
        FilterMain(argc - 1, &(argv[1]));
    }
    else if (strcmp(argv[1], "mask") == 0)
    {
        MaskMain(argc - 1, &(argv[1]));
    }
    else if (strcmp(argv[1], "merge") == 0)
    {
        MergeMain(argc - 1, &(argv[1]));
    }
    else if (strcmp(argv[1], "rank") == 0)
    {
        RankMain(argc - 1, &(argv[1]));
    }
    else if (strcmp(argv[1], "estimate") == 0)
    {
        EstimateMain(argc - 1, &(argv[1]));
    }
    else
    {
        PrintHelper();
        throw std::invalid_argument("unknown command: " + std::string(argv[1]));
    }
    return EXIT_SUCCESS;
}
