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

sbcmd="sbatch --cpus-per-task={threads}  --mem-per-cpu={cluster.mem-per-cpu-mb}"
sbcmd+=" --nodes={cluster.nodes} --ntasks={cluster.ntasks}"
sbcmd+=" --time={cluster.time} --partition={cluster.partition}"
sbcmd+=" --out={cluster.output}"

#snakemake --unlock
#snakemake -pr --keep-going \
#             --jobs 10 --cluster-config cluster.json --cluster "$sbcmd" \
#             --latency-wait 120 all1S

cp ~/atacseq-snakemake/snakemake_sub.sh /work/t/tekeller/atac_toxo
cp ~/atacseq-snakemake/cluster.json /work/t/tekeller/atac_toxo
cp ~/atacseq-snakemake/config.yml /work/t/tekeller/atac_toxo
cp ~/atacseq-snakemake/Snakefile /work/t/tekeller/atac_toxo
cp ~/atacseq-snakemake/samples.json /work/t/tekeller/atac_toxo

cd /work/t/tekeller/atac_toxo

source activate snakemake
snakemake --unlock
snakemake -j8 -pr --cluster "sbatch --time={cluster.time} --partition={cluster.partition} --mem={cluster.mem} --cpus-per-task={cluster.cpus} --out={cluster.output}" --cluster-config cluster.json --keep-going --latency-wait 120 --rerun-incomplete
