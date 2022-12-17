rm(list = ls())

library(Biostrings)
library(polyester)
library(stringr)

cmdArgs <- commandArgs(trailingOnly = TRUE)

run_seed <- 91400
fasta.path <- cmdArgs[1]
# fasta.path <- "working-bench/workdir/gencode.v34.transcripts.sub.fa"
out.dir <- cmdArgs[2]
# out.dir <- "working-bench/workdir/output/polyester_res/"
mean_read_depth <- cmdArgs[3]
# mean_read_depth <- 0.1

ref <- readDNAStringSet(fasta.path)
num_tx <- length(ref)
fold_change <- matrix(c(rep(1, num_tx), runif(n = num_tx, min = 0, max = 2)), 
                      ncol = 2, nrow = num_tx) # fold change between 0-2, uniform distribution
read_per_tx <- round(mean_read_depth * width(ref) / 100) # number of read per transcript
read_per_tx[read_per_tx == 0] <- 1 # set minimum read per transcript as 1

res.dir <- paste0(out.dir, "depth_", mean_read_depth)
dir.create(res.dir, recursive = T)
simulate_experiment(fasta = fasta.path, outdir = res.dir,
                    reads_per_transcript = read_per_tx, fold_changes = fold_change,
                    paired = TRUE, readlen = 100, num_reps = c(10, 10),
                    error_model = "uniform",
                    error_rate = 0, 
                    seed = run_seed) # error-free simulation
