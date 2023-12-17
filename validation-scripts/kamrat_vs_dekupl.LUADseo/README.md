# Run `KaMRaT` and `DE-kupl` with 130M k-mers and 50 samples

Here we test `DE-kupl`'s performance to compare to `KaMRaT`. The two jobs were running on the same CPU model indicated below.

```text
Intel(R) Xeon(R) Gold 6242 CPU @ 2.80GHz
```

## Run with `DE-kupl`

The analysis was done only including the t-test step in DE-kupl to make fair comparison.

The `Singularity` image of `DE-kupl` was built by the following command:

```bash
singularity build dekupl-run.simg docker://transipedia/dekupl-run:1.3.5
```

Then the job was launched by

```bash
DEKUPL="dekupl-run.simg"
TIME='/usr/bin/time -v'

echo "Running DE-kupl..."
BIND="-B $WKDIR/seo50sup100/:/sif_data/ -B $PWD/:/sif_pwd/ -B $WKDIR/dekupl_res/:/sif_out/"
bash -c "$TIME singularity exec $BIND $DEKUPL Rscript /dekupl/bin/Ttest_diff_method.R /dekupl/bin/TtestFilter /sif_data/seo50sup100.tab.tsv.gz /sif_data/sample_conditions_full.tsv 0.05 1 normal tumor 1 1000000 /sif_out/dekupl_tmp /sif_out/diff-counts.tsv.gz /sif_out/raw_pvals.txt.gz /sif_out/log.txt fixed"
```

## Run with `KaMRaT`

```bash
KAMRAT="KaMRaT.sif"
TIME='/usr/bin/time -v'

BIND="-B $WKDIR/seo50sup100/:/sif_data/ -B $WKDIR/kamrat_res/:/sif_out/"
echo "Running KaMRaT index..."
bash -c "$TIME singularity exec $BIND $KAMRAT kamrat index -intab /sif_data/seo50sup100.tab.tsv.gz -outdir /sif_out/index -klen 31 -unstrand -nffile /sif_data/seo50sup100.NF"
echo -e "\nRunning KaMRaT score..."
bash -c "$TIME singularity exec $BIND $KAMRAT kamrat score -idxdir /sif_out/index -scoreby ttest.padj -design /sif_data/seo50sup100.samples -outpath /sif_out/scored-counts.tsv -seltop 5000000 -withcounts"
```
