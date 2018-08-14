#! /bin/bash

#SBATCH --job-name=snakemake_gen
#SBATCH --time=120:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=8G
#SBATCH -o output.%j.%N.txt
#SBATCH -e error.%j.%N.txt

# this file is snakemake_sub.sh
#submit with
# sbatch --cpus-per-task=2 --mem=8g snakemake_sub.sh
#taken from biowulf example
#https://hpc.nih.gov/apps/snakemake.html

sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
sbcmd+=" --nodes={cluster.nodes} --ntasks={cluster.ntasks}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
sbcmd+=" --out={cluster.output}"

snakemake --unlock
snakemake -pr --keep-going \
             --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
             --latency-wait 120 all
