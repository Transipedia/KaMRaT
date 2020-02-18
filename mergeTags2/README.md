# mergeTags2 READ-ME

## DE-kupl mergeTags Idea & Interest

mergeTags is part of the DE-kupl package, which aims to assemble a given list
of k-mers into longer sequences (contigs) according to their overlap.

There are two principal interests for doing mergeTags:
  1. mergeTags groups k-mers that potentially overlap the same event, and thus
  eliminates redundant k-mer features after a differential expression analysis
  or a machine learning procedure;
  2. mergeTags extends k-mers into longer contigs, and thus provides user with
  sequence context around an event for biological and medical interpretation
  (downstream annotation or alignment).

## mergeTags2 Vision

Completeness and correctness are two main goals of k-mer assembly; however,
these may be incompatible, as a longer contig is more likely to be a wrong one.

In the first version of mergeTags (used in PMID:29284518), we concentrated 
more on the completeness by using overlap-graph tolerating 15-30 overlaps 
during assembly, while the only control of correctness was to reject assembly 
when meeting ambiguities. Based on our rough evaluation, about 20% of contigs
were mis-assembled by this approach.

In mergeTags2, we would like to achieve a better balance between completeness
and correctness by taking into consideration k-mer counts across samples
during assembly.

Moreover, mergeTags2 provides two supplementary functions:
  1. option to use the mean counts of constituent k-mers or the counts of the 
  representative k-mer (highest DE P-value) as the counts of a contig;
  2. can be applied either before or after differential analysis 
  (raw-counts or masked-counts that are without 4 columns of pvalue, meanA, 
  meanB, log2FC).

## mergeTags2 Manual

    Usage:   mergeTags2 [options] <counts.tsv>

    Options: -k INT      length of k-mers (max_value: 32) [31]  
             -m INT      min assembly overlap (max_value: k) [15]  
             -n          non-stranded merging procedure  
             -i STRING   intervention (none, pearson, spearman, contrast) [none]  
             -q STRING   quantification mode (rep, mean) [rep]  
             -s          skip differential result columns [false]

## Q&As

Question 1: What is the principle of using k-mer counts to detect mis-assemblies ?

Our basic hypothesis is that homologous k-mers should have related counts among 
samples. We propose for now 3 methods to evaluate the relationship between k-mers' 
counts among samples:

- Pearson correlation
- Spearman correlation
- Mean value of absolute contrast

Question 2: How did we evaluate the method ?

Experiment: for a given list of pre-labeled correct or incorrect contigs
created via simulation, we evaluated the performance of Pearson correlation,
Spearman correlation, and mean contrast by means of re-identifying contigs'
correctness with their quantification of the k-mers at extremities.  

Finding: The AUCs are all above 0.85, and can be as good as 0.98.  
Conclusion: Quantification is informative for finding mis-assemblies. The
performance of methods are: pearson > spearman > contrast (still very good).  

Question 3: Do the interventions have a stable performance across datasets ?

Experiment: We estimate the best thresholds of each method on 4 different
datasets.  

Finding & Conclusion: contrast has the most stable best threshold among the
methods across different datasets.

Question 4: Do the results of mergeTags2 from different intervention methods
have a good coherence/reproducibility ?

Experiment: We compare mergeTags2 results from a same differentially expressed
k-mer table and through different intervention methods. (ongoing)

Question 5: Do the calculated quantification really reveal the reality ?

Experiment: Compare the results with Reindeer. (ongoing)
