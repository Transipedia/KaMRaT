#include <zlib.h>
#include <inttypes.h>
#include <math.h> // pow()

#include <cstdlib> // free
#include <cstdint>
#include <string>

#include "kstring.h"
#include "kseq.h"
#include "khash.h"
#include "kvec.h"
#include "dna.h"
#include "stats.h"

#define VERSION "0.0.3"
#define NB_SAMPLES (nb_cols - n_diff_cols)

#define MIN_ASSEMBLY_K 15

KSEQ_INIT(gzFile, gzread)

typedef struct
{
    uint32_t ref_kmer_id;
    int nb_merged_kmers;
    char *seq;
    uint16_t ref_pos;
    float pvalue;
    float meanA, meanB;
    float log2fc;
    float *counts;
} assembly_t;

typedef struct
{
    uint32_t assembly_id;
    uint8_t revcomp;
} assembly_kmer_t;

void assembly_destroy(assembly_t *a)
{
    if (a->seq)
        free(a->seq);
    if (a->counts)
        free(a->counts);
    free(a);
}

int cmp_assembly(const void *a, const void *b)
{
    const assembly_t *ass_a = *(const assembly_t **)a;
    const assembly_t *ass_b = *(const assembly_t **)b;
    if (ass_a->ref_kmer_id == ass_b->ref_kmer_id)
    {
        return 0;
    }
    else
    {
        return (ass_a->ref_kmer_id - ass_b->ref_kmer_id);
    }
}

// Init k-mers hash
KHASH_MAP_INIT_INT64(kmers, assembly_kmer_t)

typedef khash_t(kmers) kmers_hash_t;
typedef kvec_t(assembly_t *) assemblies_array_t;

void add_assembly_kmer(kmers_hash_t *h, uint64_t kmer, uint32_t assembly_id, uint8_t rc)
{
    khiter_t k = kh_get(kmers, h, kmer);
    int dret = 0;
    if (k == kh_end(h))
    {
        k = kh_put(kmers, h, kmer, &dret);
        kh_value(h, k).assembly_id = assembly_id + 1;
        kh_value(h, k).revcomp = rc;
    }
    else
    {
        kh_value(h, k).assembly_id = 0;
    }
}

assemblies_array_t *assemble_kmers(assemblies_array_t *assembly_array,
                                   const uint8_t k_length,
                                   const uint8_t min_assembly_k,
                                   const bool stranded,
                                   const size_t n_sample,
                                   const std::string &intervention,
                                   const std::string &quant_mode,
                                   const bool rep_by_order)
{
    khiter_t k, k2;
    assemblies_array_t *a = assembly_array;

    for (uint8_t i = k_length - 1; i >= min_assembly_k; i--)
    {
        int has_merge_assemblies = 1;
        while (has_merge_assemblies)
        {
            has_merge_assemblies = 0;

            khash_t(kmers) *left_h = kh_init(kmers);
            khash_t(kmers) *right_h = kh_init(kmers);
            assemblies_array_t *new_a = (assemblies_array_t *)calloc(1, sizeof(assemblies_array_t));
            kv_init(*new_a);

            // fprintf(stderr, "Merging sequences with k = %d\n", i);

            // fprintf(stderr, "Creating kmer-end indexes\n");

            // First we index the start and ends of each assembly
            for (size_t j = 0; j < kv_size(*a); j++)
            {
                assembly_t *assembly = kv_A(*a, j);
                if (strlen(assembly->seq) <= i)
                    continue;
                uint64_t left_kmer = dna_to_int(assembly->seq, i);
                uint64_t right_kmer = dna_to_int(&assembly->seq[strlen(assembly->seq) - i], i);
                uint8_t is_left_kmer_rc = 0, is_right_kmer_rc = 0;

                // If we are not in stranded mode, we index the canonical k-i k-mers
                if (!stranded)
                {
                    uint64_t left_kmer_rc = int_revcomp(left_kmer, i);
                    uint64_t right_kmer_rc = int_revcomp(right_kmer, i);
                    if (left_kmer_rc < left_kmer)
                    {
                        left_kmer = left_kmer_rc;
                        is_left_kmer_rc = 1;
                    }
                    if (right_kmer_rc < right_kmer)
                    {
                        right_kmer = right_kmer_rc;
                        is_right_kmer_rc = 1;
                    }
                }

                if (is_left_kmer_rc)
                {
                    add_assembly_kmer(right_h, left_kmer, j, is_left_kmer_rc);
                }
                else
                {
                    add_assembly_kmer(left_h, left_kmer, j, is_left_kmer_rc);
                }

                if (is_right_kmer_rc)
                {
                    add_assembly_kmer(left_h, right_kmer, j, is_right_kmer_rc);
                }
                else
                {
                    add_assembly_kmer(right_h, right_kmer, j, is_right_kmer_rc);
                }
            }

            // fprintf(stderr, "Merging sequences\n");
            // Now we try to merge assemblies
            for (k = kh_begin(right_h); k != kh_end(right_h); ++k)
            {
                if (!kh_exist(right_h, k))
                    continue;

                if (kh_value(right_h, k).assembly_id == 0)
                    continue;

                assembly_t *left_assembly = kv_A(*a, kh_value(right_h, k).assembly_id - 1);

                if (!left_assembly)
                    continue;

                k2 = kh_get(kmers, left_h, kh_key(right_h, k));

                if (k2 != kh_end(left_h) && kh_value(left_h, k2).assembly_id != 0 && kh_value(left_h, k2).assembly_id != kh_value(right_h, k).assembly_id)
                {
                    // Merge the two sequences into a new assembly
                    assembly_t *right_assembly = kv_A(*a, kh_value(left_h, k2).assembly_id - 1);

                    if (!right_assembly)
                        continue;

                    float pearson = calc_pearson_corr(left_assembly->counts, right_assembly->counts, n_sample),
                          spearman = calc_spearman_corr(left_assembly->counts, right_assembly->counts, n_sample),
                          mean_ctrst = calc_mean_contrast(left_assembly->counts, right_assembly->counts, n_sample);

                    if (intervention == MET_PEARSON && pearson < MIN_PEARSON)
                    {
                        continue;
                    }
                    else if (intervention == MET_SPEARMAN && spearman < MIN_SPEARMAN)
                    {
                        continue;
                    }
                    else if (intervention == MET_CONTRAST && mean_ctrst > MAX_CONTRAST)
                    {
                        continue;
                    }

                    has_merge_assemblies = 1;

                    assembly_t *assemblies_merge, *removed_assembly;
                    // Keep the k-mers with lowest pvalue to represent the assembly
                    if (rep_by_order)
                    {
                        if (left_assembly->ref_kmer_id < right_assembly->ref_kmer_id)
                        {
                            assemblies_merge = left_assembly;
                            removed_assembly = right_assembly;
                        }
                        else
                        {
                            assemblies_merge = right_assembly;
                            removed_assembly = left_assembly;
                            assemblies_merge->ref_pos += (strlen(left_assembly->seq) - i);
                        }
                    }
                    else
                    {
                        if (left_assembly->pvalue < right_assembly->pvalue)
                        {
                            assemblies_merge = left_assembly;
                            removed_assembly = right_assembly;
                        }
                        else
                        {
                            assemblies_merge = right_assembly;
                            removed_assembly = left_assembly;
                            assemblies_merge->ref_pos += (strlen(left_assembly->seq) - i);
                        }
                    }

                    // Create a new string for the assembly merging
                    char *merge_seq = (char *)malloc(strlen(right_assembly->seq) + strlen(left_assembly->seq) - i + 1);

                    // If the assemblies are not in the same orientation, we reverse one of them
                    if (kh_value(left_h, k2).revcomp != kh_value(right_h, k).revcomp)
                    {
                        if (kh_value(left_h, k2).revcomp)
                        {
                            char *buffer = (char *)malloc(strlen(right_assembly->seq) + 1);
                            revcomp(right_assembly->seq, buffer, strlen(right_assembly->seq));
                            free(right_assembly->seq);
                            right_assembly->seq = buffer;
                        }
                        else
                        {
                            char *buffer = (char *)malloc(strlen(left_assembly->seq) + 1);
                            revcomp(left_assembly->seq, buffer, strlen(left_assembly->seq));
                            free(left_assembly->seq);
                            left_assembly->seq = buffer;
                        }
                    }
                    else if (kh_value(left_h, k2).revcomp && kh_value(right_h, k).revcomp)
                    {
                        // If both assemblies are in RC, we swap them
                        assembly_t *tmp_assembly = left_assembly;
                        left_assembly = right_assembly;
                        right_assembly = tmp_assembly;
                    }

                    strcpy(merge_seq, left_assembly->seq);
                    strcat(merge_seq, &right_assembly->seq[i]);

                    // std::cerr << merge_seq << "\t" << left_assembly->seq << "\t" << right_assembly->seq;

                    free(right_assembly->seq);
                    free(left_assembly->seq);
                    right_assembly->seq = NULL;
                    left_assembly->seq = NULL;
                    assemblies_merge->seq = merge_seq;
                    if (quant_mode == "mean")
                    {
                        for (size_t s = 0; s < n_sample; s++)
                        {
                            float left_count_s = left_assembly->counts[s] * left_assembly->nb_merged_kmers,
                                  right_count_s = right_assembly->counts[s] * right_assembly->nb_merged_kmers;
                            assemblies_merge->counts[s] = (left_count_s + right_count_s) / (left_assembly->nb_merged_kmers + right_assembly->nb_merged_kmers);
                            // std::cerr << "\t" << assemblies_merge->counts[s];
                        }
                        // std::cerr << std::endl;
                    }
                    assemblies_merge->nb_merged_kmers = left_assembly->nb_merged_kmers + right_assembly->nb_merged_kmers;
                    free(removed_assembly->counts);
                    removed_assembly->counts = NULL;
                    assembly_destroy(removed_assembly);

                    kv_push(assembly_t *, *new_a, assemblies_merge);
                    kv_A(*a, kh_value(right_h, k).assembly_id - 1) = NULL;
                    kv_A(*a, kh_value(left_h, k2).assembly_id - 1) = NULL;
                }
            }

            // fprintf(stderr, "Add un-merged sequences\n");
            // Add assemblies that were not merged
            for (size_t j = 0; j < kv_size(*a); j++)
            {
                assembly_t *assembly = kv_A(*a, j);
                if (assembly)
                {
                    kv_push(assembly_t *, *new_a, assembly);
                }
            }
            // fprintf(stderr, "%zu assemblies after merging\n", kv_size(*new_a));
            // Delete previous assebly array and replace it with the new one
            kv_destroy(*a);
            a = new_a;
            // Delete hashes
            kh_destroy(kmers, left_h);
            kh_destroy(kmers, right_h);
        }
    }
    return a;
}

int main(int argc, char *argv[])
{
    char *counts_file;
    uint8_t k_length = 31;
    bool stranded = true, rep_by_order(false);
    std::string intervention = "none", quant_mode = "rep";
    uint8_t min_assembly_k = MIN_ASSEMBLY_K;
    uint8_t n_diff_cols = 4;

    int c;
    while ((c = getopt(argc, argv, "nk:m:i:q:s:p")) >= 0)
    {
        switch (c)
        {
        case 'n':
            stranded = false;
            break;
        case 'k':
            k_length = atoi(optarg);
            break;
        case 'm':
            min_assembly_k = atoi(optarg);
            break;
        case 'i':
            intervention = optarg;
            break;
        case 'q':
            quant_mode = optarg;
            break;
        case 's':
            n_diff_cols = 0;
            break;
        case 'p':
            rep_by_order = true;
            break;
        }
    }

    if (optind == argc)
    {
        std::cerr << std::endl;
        std::cerr << "Usage:   mergeTags2 [options] <counts.tsv>" << std::endl
                  << std::endl;
        std::cerr << "Options: -k INT        length of k-mers (max_value: 32) [" << int(k_length) << "]" << std::endl;
        std::cerr << "         -m INT        min assembly overlap (max_value: k) [" << int(min_assembly_k) << "]" << std::endl;
        std::cerr << "         -n            unstranded merging procedure" << std::endl;
        std::cerr << "         -i STRING     intervention method ("
                  << MET_NONE << ", " << MET_PEARSON << ", " << MET_SPEARMAN << ", " << MET_CONTRAST << ") [" << MET_NONE << "]" << std::endl;
        std::cerr << "         -q STRING     quantification mode (rep, mean) [rep]" << std::endl;
        std::cerr << "         -s            skip differential result columns [false]" << std::endl;
        std::cerr << "         -p BOOLEAN    representative k-mer by input order [false]" << std::endl
                  << std::endl;
        return EXIT_SUCCESS;
    }

    counts_file = argv[optind++];

    std::cerr << "========> mergeTags-2 RUN INFO <========" << std::endl;
    std::cerr << "Stranded: " << (stranded ? "true" : "false") << std::endl;
    std::cerr << "Intervention method: " << intervention;
    if (intervention == MET_PEARSON)
    {
        std::cerr << " [threshold = " << MIN_PEARSON << "]" << std::endl;
    }
    else if (intervention == MET_SPEARMAN)
    {
        std::cerr << " [threshold = " << MIN_SPEARMAN << "]" << std::endl;
    }
    else if (intervention == MET_CONTRAST)
    {
        std::cerr << " [threshold = " << MAX_CONTRAST << "]" << std::endl;
    }
    else
    {
        std::cerr << std::endl;
    }
    std::cerr << "Quantification mode: " << quant_mode << std::endl;

    // fprintf(stderr, "Loading k-mers into memory\n");

    size_t nb_cols = 0;
    int nb_kmers = 0, dret = 0;
    gzFile fp;
    kstream_t *ks;
    kstring_t *str;
    assemblies_array_t *up_assemblies = (assemblies_array_t *)calloc(1, sizeof(assemblies_array_t));
    assemblies_array_t *down_assemblies = (assemblies_array_t *)calloc(1, sizeof(assemblies_array_t));

    kv_init(*up_assemblies);
    kv_init(*down_assemblies);
    str = (kstring_t *)calloc(1, sizeof(kstring_t));

    fp = gzopen(counts_file, "r");
    if (!fp)
    {
        // fprintf(stderr, "Failed to open %s\n", counts_file);
        exit(EXIT_FAILURE);
    }

    ks = ks_init(fp);
    // Read and Print Header Line
    std::cout << "nb_merged_kmers\tcontig";
    while (ks_getuntil(ks, 0, str, &dret) > 0 && dret != '\n')
    {
        std::cout << "\t" << str->s;
        nb_cols++;
    }
    std::cout << "\t" << str->s << std::endl;

    std::cerr << "Total sample number: " << NB_SAMPLES << std::endl;

    while (ks_getuntil(ks, 0, str, &dret) >= 0)
    {
        assembly_t *assembly = (assembly_t *)calloc(1, sizeof(assembly_t));
        assembly->ref_kmer_id = nb_kmers;
        assembly->seq = ks_release(str);
        assembly->ref_pos = 0;
        assembly->nb_merged_kmers = 1;
        assembly->counts = (float *)calloc(NB_SAMPLES, sizeof(float));

        // Get all remaining columns
        for (size_t i = 0; i < nb_cols; i++)
        {
            if (ks_getuntil(ks, 0, str, &dret) > 0)
            {
                if (i < n_diff_cols && i == 0)
                {
                    assembly->pvalue = atof(str->s);
                }
                else if (i < n_diff_cols && i == 1)
                {
                    assembly->meanA = atof(str->s);
                }
                else if (i < n_diff_cols && i == 2)
                {
                    assembly->meanB = atof(str->s);
                }
                else if (i < n_diff_cols && i == 3)
                {
                    assembly->log2fc = atof(str->s);
                }
                else if (i >= n_diff_cols)
                {
                    assembly->counts[i - n_diff_cols] = atof(str->s);
                }
                else
                {
                    std::cerr << "An unconsidered condition !" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                // fprintf(stderr, "Invalid number of columns at line %d\n", nb_kmers + 2);
                exit(EXIT_FAILURE);
            }
        }

        if (dret == '\n')
        {
            if (assembly->log2fc > 0)
            {
                kv_push(assembly_t *, *up_assemblies, assembly);
            }
            else
            {
                kv_push(assembly_t *, *down_assemblies, assembly);
            }
        }
        else
        {
            fprintf(stderr, "Invalid number of columns at line %d\n", nb_kmers + 2);
            exit(EXIT_FAILURE);
        }

        nb_kmers++;
    }
    ks_destroy(ks);
    gzclose(fp);
    free(str->s);
    free(str);

    // fprintf(stderr, "%d k-mers loaded\n", nb_kmers);
    // fprintf(stderr, "Assembling k-mers\n");

    // Assemble kmers with positive log2fc
    // fprintf(stderr, "Assembling %d up regulated k-mers\n", kv_size(*up_assemblies));
    up_assemblies = assemble_kmers(up_assemblies, k_length, min_assembly_k, stranded, NB_SAMPLES, intervention, quant_mode, rep_by_order);
    // Assemble kmers with negative log2fc
    // fprintf(stderr, "Assembling %d down regulated k-mers\n", kv_size(*down_assemblies));
    down_assemblies = assemble_kmers(down_assemblies, k_length, min_assembly_k, stranded, NB_SAMPLES, intervention, quant_mode, rep_by_order);
    // Merge all assemblies together
    // fprintf(stderr, "Merging up and down assemblies\n");
    for (size_t j = 0; j < kv_size(*down_assemblies); j++)
    {
        kv_push(assembly_t *, *up_assemblies, kv_A(*down_assemblies, j));
    }
    kv_destroy(*down_assemblies);

    // Sort assemblies by ref_kmer_id
    qsort(up_assemblies->a, kv_size(*up_assemblies), sizeof(up_assemblies->a[0]), cmp_assembly);

    // Print assembled contigs
    for (size_t i = 0; i < kv_size(*up_assemblies); i++)
    {
        assembly_t *assembly = kv_A(*up_assemblies, i);
        std::string tmp_seq(assembly->seq);
        std::cout << assembly->nb_merged_kmers << "\t" << assembly->seq << "\t" << tmp_seq.substr(assembly->ref_pos, k_length);
        if (n_diff_cols > 0)
        {
            std::cout << "\t";
            std::cout.precision(6);
            std::cout << std::fixed
                      << assembly->pvalue << "\t"
                      << assembly->meanA << "\t"
                      << assembly->meanB << "\t"
                      << assembly->log2fc;
        }
        std::cout.precision(2);
        std::cout << std::fixed;
        for (size_t j = 0; j < NB_SAMPLES; j++)
        {
            std::cout << "\t" << assembly->counts[j];
        }
        std::cout << std::endl;
    }
    kv_destroy(*up_assemblies);

    //free(a);
    return 0;
}
