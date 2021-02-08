rm(list = ls())

library(Biostrings)
library(polyester)
library(stringr)
library(magrittr)
library(ggplot2)

cmdArgs <- commandArgs(trailingOnly = TRUE)
fasta.path <- cmdArgs[1]
out.dir <- cmdArgs[2]

# fasta.path <- "/home/haoliang.xue/media/data/kamrat/paper/c_simulated_reads/a_ref/ref4simu.fa"
# out.dir <- "/home/haoliang.xue/media/data/kamrat/paper/c_simulated_reads/d_witherr_100vs100"

nb.rep <- 100
run_seed <- 91400

set.seed(run_seed)

ref.fasta <- readDNAStringSet(fasta.path)

## For Transcripts
tx.fasta <- ref.fasta[str_detect(names(ref.fasta), pattern = "ENST")]
tx.nb <- length(tx.fasta) # 3765 transcripts
tx.nread <- round(20 * width(tx.fasta) / 100) # ~20x coverage for each transcript

tx.fc <- matrix(c(rep(1, tx.nb),
                  runif(tx.nb, 
                        min = ifelse(tx.nread < 30, yes = 3, no = 0.2),
                        max = ifelse(tx.nread < 30, yes = 5, no = 1.8))),
                ncol = 2,
                byrow = FALSE)

tx.mean <- tx.nread * tx.fc     # mu parameter in NB model
tx.size <- tx.nread * tx.fc / 3 # size parameter in NB model = read_number * fold_change / 3

tx.count.mat <- lapply(1 : tx.nb,
                      function(i) return(c(rnbinom(n = nb.rep, size = tx.size[i, 1], mu = tx.mean[i, 1]),
                                           rnbinom(n = nb.rep, size = tx.size[i, 2], mu = tx.mean[i, 2])))) %>%
    do.call(what = rbind) # negative binomial with different means between conditions

grp1.means <- rowMeans(tx.count.mat[, 1 : nb.rep])
grp2.means <- rowMeans(tx.count.mat[, (nb.rep + 1) : (2 * nb.rep)])
simu.fc <- grp2.means / grp1.means

rownames(tx.count.mat) <- names(tx.fasta)

ggplot() + 
    geom_histogram(aes(x = simu.fc), binwidth = 0.02, fill = "navy", alpha = 0.5) +
    xlab("fold-change value") +
    theme(text = element_text(size = 22)) +
    ggsave(paste0(out.dir, "/transcripts_fc_histo.png"), width = 21, height = 9)
rm(tx.fc, tx.mean, tx.size, grp1.means, grp2.means, simu.fc, tx.nb, tx.nread)

## For Variation Events
ev.fasta <- ref.fasta[str_detect(names(ref.fasta), pattern = "rs")]
ev.nb <- length(ev.fasta)
ev.nread <- round(3 * nchar(ev.fasta) / 100) # ~3x coverage for variation events

ev.pres_grp <- sample(c(1, 2), size = ev.nb, replace = T)
ev.mean <- ev.nread
ev.size <- ev.nread / 3

ev.count.mat <- lapply(1 : ev.nb,
                       function (i) {
                           count.row <- rep(0, 2 * nb.rep)
                           express.smp.nb <- round(runif(1, min = 10, max = 20))
                           express.smp.i <- sample(1 : nb.rep, size = express.smp.nb, replace = F)
                           count.row[express.smp.i] <- rnbinom(express.smp.nb, size = ev.size[i], mu = ev.mean[i])
                           return(count.row) 
                       }) %>%
    do.call(what = rbind)

rownames(ev.count.mat) <- names(ev.fasta)

rm(ev.pres_grp, ev.nb, ev.nread, ev.mean, ev.size)

## Combine Count Matrices and Simulate
simu.fasta <- c(tx.fasta, ev.fasta)
simu.fasta.path <- paste0(out.dir, "/ref4simu.fa")
writeXStringSet(filepath = simu.fasta.path, x = simu.fasta)
simu.count.mat <- rbind(tx.count.mat, ev.count.mat)
write.table(simu.count.mat, paste0(out.dir, "/countmat4simu.tsv"), row.names = T, col.names = T, quote = F, sep = "\t")
simulate_experiment_countmat(fasta = simu.fasta.path, readmat = simu.count.mat, outdir = out.dir,
                             paired = TRUE, readlen = 100, error_model = "illumina4",
                             seed = run_seed)
