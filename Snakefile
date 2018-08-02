from os.path import join
GENOME='~/toxo_kim/genomes/hg38.fa'
fastq_dirs='/work/t/tekeller/atac_toxo'


FILES = json.load(open(config['SAMPLES_JSON']))

import csv
import os

SAMPLES = sorted(FILES.keys())


rule all:
    input:
        expand("{sample}.txt",sample=SAMPLES)

rule clean_fastq:
    input:
        genome=GENOME,
        fwd=expand("/work/t/tekeller/atac_toxo/{sample}.1_val_1.fq.gz",run=RUN,sample=SAMPLES),
        rev=expand("/work/t/tekeller/atac_toxo/{sample}.2_val_2.fq.gz",run=RUN,sample=SAMPLES)
    output:
       expand("{sample}.txt",sample=SAMPLES)
    shell:
        """
        echo {input.genome} {input.fwd} {input.rev} > {output}
        """

