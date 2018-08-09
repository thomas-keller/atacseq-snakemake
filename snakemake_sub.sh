#! /bin/bash
# this file is snakemake_sub.sh
#submit with
# sbatch --cpus-per-task=2 --mem=8g snakemake.sh
#taken from biowulf example
#https://hpc.nih.gov/apps/snakemake.html

sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
sbcmd+=" --out={cluster.output}"
snakemake -pr --keep-going --local-cores $SLURM_CPUS_PER_TASK \
             --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
             --latency-wait 120 all

snakemake -pr --keep-going --local-cores $SLURM_CPUS_PER_TASK \
    --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
    --latency-wait 120 all