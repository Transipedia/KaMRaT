# Evaluate `DE-kupl` Performance with 130M k-mers and 50 samples

Here we test `DE-kupl`'s performance to compare to `KaMRaT`.

The analysis was done inside the `DE-kupl`'s official singularity image.

```bash
singularity build dekupl-run.simg docker://transipedia/dekupl-run:1.3.5
```

Then the `DE-kupl` analysis was run with the given config file and modified Snakefile using following commands:

```bash
singularity shell -c -W $PWD -B seo50sup100/:/sif_data/ -B /usr/bin/:/bin4time/ -B $PWD/:/sif_pwd/ -B dekupl_res/:/sif_out/ dekupl-run.simg
cd /sif_pwd
bash -c '/bin4time/time -v snakemake --configfile /sif_pwd/my-config.json -s /sif_pwd/Snakefile --cores 8' 2>&1 | tee /sif_pwd/dekupl.logerr
```

