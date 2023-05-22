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
mean_read_depth <- as.numeric(cmdArgs[3])
# mean_read_depth <- 0.1
mode <- cmdArgs[4]

cat(fasta.path, out.dir, mean_read_depth, mode, sep = "\n")

ref <- readDNAStringSet(fasta.path)
num_tx <- length(ref)
read_per_tx <- round(mean_read_depth * width(ref) / 100) # number of read per transcript
read_per_tx[read_per_tx == 0] <- 1 # set minimum read per transcript as 1
if (mode == "err-free") {
	res.dir <- paste0(out.dir, "err-free/depth_", mean_read_depth, "/", sub(".fa", "", basename(fasta.path)))
	simulate_experiment(fasta = fasta.path, outdir = res.dir,
		    	    reads_per_transcript = read_per_tx,
			    fold_changes = rep(1, num_tx), num_reps = c(1, 1),
    			    paired = TRUE, readlen = 100,
	    		    error_model = "uniform", error_rate = 0,
			    seed = run_seed)
} else {
	res.dir <- paste0(out.dir, "err-illumina5/depth_", mean_read_depth, "/", sub(".fa", "", basename(fasta.path)))
	simulate_experiment(fasta = fasta.path, outdir = res.dir,
			    reads_per_transcript = read_per_tx,
			    fold_changes = rep(1, num_tx), num_reps = c(1, 1),
			    paired = TRUE, readlen = 100,
			    error_model = "illumina5",
			    seed = run_seed)
}
cat("Polyester completed!\n")
