rm(list = ls())

library(Biostrings)
library(polyester)
library(stringr)

cmdArgs <- commandArgs(trailingOnly = TRUE)
fasta.path <- cmdArgs[1]
out.dir <- cmdArgs[2]
with.err <- cmdArgs[3]

run_seed <- 91400

# fasta.path <- "/home/haoliang.xue/media/data/kamrat/paper/b_reference_gen/a_selected_transcripts.fa"
# out.dir <- "/home/haoliang.xue/media/data/kamrat/paper/c_simulated_reads/error_free"
# with.err <- false

num_tx <- count_transcripts(fasta.path)

ref <- readDNAStringSet(fasta.path)

fold_change <- matrix(1, ncol = 2, nrow = num_tx)
read_per_tx <- round(20 * width(ref) / 100) # ~20x coverage for each transcript

# remove quotes from transcript IDs:
if (with.err) {	
	simulate_experiment(fasta = fasta.path, 
        	            outdir = out.dir,
                	    num_reps = c(5, 5),
	                    reads_per_transcript = read_per_tx,
        	            fold_changes = fold_change,
                	    paired = TRUE,
	                    readlen = 100,
        	            error_model = "illumina4",
	                    seed = run_seed)
} else {
	simulate_experiment(fasta = fasta.path, 
        	            outdir = out.dir,
                	    num_reps = c(5, 5),
	                    reads_per_transcript = read_per_tx,
        	            fold_changes = fold_change,
                	    paired = TRUE,
	                    readlen = 100,
        	            error_model = "uniform",
                	    error_rate = 0, 
	                    seed = run_seed)
}
