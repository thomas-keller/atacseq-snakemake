#! /bin/bash
# this file is snakemake_sub.sh
#submit with
# sbatch --ntasks=2 --mem=8g snakemake_sub.sh
#taken from biowulf example
#https://hpc.nih.gov/apps/snakemake.html

sbcmd="sbatch --ntasks={threads} --mem={cluster.mem}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
sbcmd+=" --out={cluster.output}"
snakemake -pr --unlock --keep-going --local-cores $SLURM_TASKS_PER_NODE \
             --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
             --latency-wait 120 all

snakemake -pr --keep-going --local-cores $SLURM_TASKS_PER_NODE \
    --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
    --latency-wait 120 all