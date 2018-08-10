#!/bin/bash

#SBATCH --job-name=snakemake_gen
#SBATCH --time=08:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=8G
#SBATCH -o output.%j.%N.txt
#SBATCH -e error.%j.%N.txt

#
# check for clean
#
# https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-do-i-remove-all-files-created-by-snakemake-i-e-like-make-clean
if [ "$1" == "clean" ]; then
    echo 'rm $(snakemake --summary | tail -n+2 | cut -f1)'
    snakemake --summary | tail -n+2 | cut -f1
    rm -f $(snakemake --summary | tail -n+2 | cut -f1)
    exit 0
fi

#
# launch snakemake to run jobs on UAB CHEAHA via SLURM
#
SM_PARAMS="job-name ntasks partition time mail-user mail-type error output"
SM_ARGS="--cpus-per-task {threads} --mem {cluster.mem}"
for P in ${SM_PARAMS}; do SM_ARGS="$SM_ARGS --$P={cluster.$P}"; done
echo "SM_ARGS: ${SM_ARGS}"

# our SLURM error/output paths expect a logs/ subdir in PWD
mkdir -p logs

snakemake \
    $* \
     --latency-wait 30 \
    -j 999 \
    --cluster-config cluster.json \
--cluster "sbatch $SM_ARGS"