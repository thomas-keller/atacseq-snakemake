from os.path import join
GENOME='~/toxo_kim/hg38.fa'
fastq_dirs='/work/t/tekeller/atac_toxo'
RUN,SAMPLES=glob_wildcards("/work/t/tekeller/atac_toxo/{run}/{sample}.1_val_1.fq.gz")

PATTERN_R1 = '{sample}.1_val_1.fq.gz'
PATTERN_R2 = '{sample}.2_val_2.fq.gz'



rule all:
    input:
        "{sample}.txt"

rule clean_fastq:
    input:
        genome=GENOME,
        fwd=expand("/work/t/tekeller/atac_toxo/{run}/{sample}.1_val_1.fq.gz",run=RUN,sample=SAMPLES),
        rev=expand("/work/t/tekeller/atac_toxo/{run}/{sample}.2_val_2.fq.gz",run=RUN,sample=SAMPLES)
    output:
       "{sample}.txt"
    shell:
        """
        echo {input.genome} {input.fwd} {input.rev} > {output}
        """

