rm(list = ls())

library(Biostrings)
library(polyester)
library(stringr)

cmdArgs <- commandArgs(trailingOnly = TRUE)

run_seed <- 91400
fasta.path <- cmdArgs[1]
out.dir <- cmdArgs[2]

ref <- readDNAStringSet(fasta.path)
num_tx <- length(ref)
fold_change <- matrix(c(rep(1, num_tx), runif(n = num_tx, min = 0, max = 2)), 
                      ncol = 2, nrow = num_tx) # fold change between 0-2, uniform distribution
read_per_tx <- round(10 * width(ref) / 100) # ~10x coverage for each transcript, minimum coverage is 1

dir.create(out.dir, recursive = T)
simulate_experiment(fasta = fasta.path, outdir = out.dir,
                    reads_per_transcript = read_per_tx, fold_changes = fold_change,
                    paired = TRUE, readlen = 100, num_reps = c(10, 10),
                    error_model = "uniform",
                    error_rate = 0, 
                    seed = run_seed) # error-free simulation
