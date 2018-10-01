import json
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output
import csv

# Globals ---------------------------------------------------------------------
configfile: "config.yml"

GENOME='/home/t/tekeller/toxo_kim/genomes/hg38.fa'
TOXO='/home/t/tekeller/toxo_kim/genomes/ToxoDB-38_TgondiiME49_Genome.fasta'
fastq_dirs='/work/t/tekeller/atac_toxo'


FILES = json.load(open(config['SAMPLES_JSON']))



SAMPLES = sorted(FILES.keys())

## list all samples by sample_name and sample_type
#MARK_SAMPLES = []
#for sample in SAMPLES:
#    for sample_type in FILES[sample].keys():
#        MARK_SAMPLES.append(sample + ";" + sample_type)

CONTROL = config["control"]
CONTROLS = [sample for sample in SAMPLES if CONTROL in sample]
CASES = [sample for sample in SAMPLES if CONTROL not in sample]


CASE_FILES=expand("/work/t/tekeller/atac_toxo/{case}.1_val_1.fq.gz",case=CASES)

logdir = os.path.join(os.getcwd(), "logs/slurm")
os.makedirs(logdir, exist_ok=True)

# Functions -------------------------------------------------------------------

def rstrip(text, suffix):
    # Remove a suffix from a string.
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

# Rules -----------------------------------------------------------------------

#CONTROL_MERGED_FASTQ = expand("02seq/{control}_merged.fq.gz", control = CONTROLS)
CASE_CLEAN_HG=expand("03cln/{case}_clean_hg38.fq",case=CASES)
CASE_CLEAN_TOXO=expand("03cln/{case}_clean_ToxoDB-38_TgondiiME49_Genome.fq",case=CASES)
ALL_SAMPLES=expand("/work/t/tekeller/atac_toxo/{sample}.{read}.fq.gz",sample=SAMPLES,read=['1_val_1','2_val_2'])

FASTQC_ALL=expand('/work/t/tekeller/atac_toxo/01fqc/{sample}.{read}.fastqc.zip',sample=SAMPLES,read=['1_val_1','2_val_2'])
# ATAQV_CTL=expand("04aln/{control}.sorted.bam.ataqv.json",control=CONTROLS)
# ATAQV_HG=expand("04aln/{case}_hg.sorted.bam.ataqv.json",case=CASES)
#ATAQV_ALL=ATAQV_CTL+ATAQV_HG

ALN_HG=expand("04aln/{case}_hg.sorted.bam",case=CASES)
ALN_TOXO=expand("04aln/{case}_toxo.sorted.bam",case=CASES)
ALN_CTL=expand("04aln/{control}.sorted.bam",control=CONTROLS)
ALN_ALL=ALN_HG+ALN_TOXO+ALN_CTL


FLAG_HG=expand("04aln/{case}_hg.sorted.bam.flagstat",case=CASES)
FLAG_TOXO=expand("04aln/{case}_toxo.sorted.bam.flagstat",case=CASES)
FLAG_CTL=expand("04aln/{control}_toxo.sorted.bam.flagstat",control=CONTROLS)
FLAG_ALL=FLAG_HG+FLAG_TOXO+FLAG_CTL

RM_CHRM_CTL=expand("06aln_exclude_chrM/{control}_exclude_chrM.sorted.bam",control=CONTROLS)
RM_CHRM_CASES=expand("06aln_exclude_chrM/{case}_hg_exclude_chrM.sorted.bam",case=CASES)
RM_CHRM_ALL=RM_CHRM_CTL+RM_CHRM_CASES

#NUCL_CASE=expand("09nucleoATAC/{case}_nucleoATAC.occpeaks.bed.gz",case=CASES)
#NUCL_CTL=expand("09nucleoATAC/{control}_nucleoATAC.occpeaks.bed.gz",control=CONTROLS)
#NUCL_ALL=NUCL_CASE+NUCL_CTL

#ATAQV=["11ATAC_qc_html"]

print(CASES)
print(CONTROLS)
rule all:
    input:
		#CONTROL_MERGED_FASTQ + CASE_CLEAN_HG + CASE_CLEAN_TOXO+ALN_ALL+FLAG_ALL+NUCL_ALL+RM_CHRM_ALL
        ALL_SAMPLES+ CASE_CLEAN_HG + CASE_CLEAN_TOXO+ALN_ALL+FLAG_ALL+ATAQV

CONTROL_FILES=expand('/work/t/tekeller/atac_toxo/{control}.1_val_1.fq.gz',control=CONTROLS)
ALL_SAMPLES=expand("/work/t/tekeller/atac_toxo/{sample}.1_val_1.fq.gz",sample=SAMPLES)
ALL_CASES=expand("/work/t/tekeller/atac_toxo/{case}.1_val_1.fq.gz",case=CASES)

rule fastqc:
    input:  "/work/t/tekeller/atac_toxo/{sample}.1_val_1.fq.gz", "/work/t/tekeller/atac_toxo/{sample}.2_val_2.fq.gz"
    output:"01fqc/{sample}.1_val_1.fastqc.zip", "01fqc/{sample}.2_val_2.fastqc.zip"
    log:    "00log/{sample}_fastqc"
    threads: 1
    params : jobname = "{sample}"
    message: "fastqc {input}: {threads}"
    shell:
        """
	# fastqc works fine on .gz file as well
        module load fastqc
        fastqc -o 02fqc -f fastq --noextract {input[0]} {input[1]} 2> {log}
        """

rule trim_adapter:
 	input: "/work/t/tekeller/atac_toxo/{sample}.1_val_1.fq.gz", "/work/t/tekeller/atac_toxo/{sample}.2_val_2.fq.gz"
 	output: "02trim/{sample}.1_val_1.trimmed.fastq.gz", "02trim/{sample}.2_val_2.trimmed.fastq.gz"
 	log: "00log/{sample}_trim_adaptor.log"
 	threads: 1
 	params : jobname = "{sample}"
 	message: "trim_adaptor {input}: {threads}"
 	shell:
 		"""
 		# python2.x
 		#source activate root
		trim_adapters {input[0]} {input[1]} 2> {log}
		mv /work/t/tekeller/atac_toxo/{wildcards.sample}.1_val_1.trimmed.fastq.gz 02trim/
		mv /work/t/tekeller/atac_toxo/{wildcards.sample}.2_val_2.trimmed.fastq.gz 02trim/
 		"""


rule clean_fastq:
    input:
        genome=GENOME,
        toxo=TOXO,
        fwd="02trim/{case}.1_val_1.trimmed.fastq.gz",
        rev="02trim/{case}.2_val_2.trimmed.fastq.gz"
    output:
       hum_cl="03cln/{case}_clean_hg38.fq",
       tox_cl="03cln/{case}_clean_ToxoDB-38_TgondiiME49_Genome.fq"
    shell:
        """
        bbsplit.sh -Xmx28g -t=8 in={input.fwd} in2={input.rev} basename=03cln/{wildcards.case}_clean_%.fq
        """

#following pyflow-ATACseq
#https://github.com/crazyhottommy/pyflow-ATACseq
## the later step will remove chrM from the bam file and coordinate sort the bam
## so I did not cooridnate sort the bam at this step to save some time.
rule align_cases_hg:
	input: "03cln/{case}_clean_hg38.fq"
	output: "04aln/{case}_hg.sorted.bam", "00log/{case}.align"
	params: jobname = "{case}"
	message: "aligning {input}: {threads} threads"
	log:
		bowtie2 = "00log/{case}_hg.align",
		markdup = "00log/{case}_hg.markdup"
	shell:
		"""
		## samblaster mark duplicates for read id grouped reads. I do not coordinate sort the bam
		bowtie2 --threads 5  -X2000 -x {config[idx_bt2]} --interleaved {input[0]} 2> {log.bowtie2} \
		| samblaster 2> {log.markdup} \
		| samtools view -Sb - > {output[0]}
		"""

rule align_cases_toxo:
	input: "/work/t/tekeller/atac_toxo/03cln/{case}_clean_ToxoDB-38_TgondiiME49_Genome.fq"
	output: "04aln/{case}_toxo.sorted.bam", "00log/{case}.align"
	params: jobname = "{case}"
	message: "aligning {input}: {threads} threads"
	log:
		bowtie2 = "00log/{case}_toxo.align",
		markdup = "00log/{case}_toxo.markdup"
	shell:
		"""
		## samblaster mark duplicates for read id grouped reads. I do not coordinate sort the bam
		bowtie2 --threads 5  -X2000 -x {config[idx_bt2_toxo]} --interleaved {input[0]} 2> {log.bowtie2} \
		| samblaster 2> {log.markdup} \
		| samtools view -Sb - > {output[0]}
		"""


rule align_control:
	input: "/work/t/tekeller/atac_toxo/{control}.1_val_1.fq.gz","/work/t/tekeller/atac_toxo/{control}.2_val_2.fq.gz"
	output: "04aln/{control}.sorted.bam", "00log/{control}.align"
	params: jobname = "{control}"
	message: "aligning {input}: {threads} threads"
	log:
		bowtie2 = "00log/{control}.align",
		markdup = "00log/{control}.markdup"
	shell:
		"""
		module add apps/samtools
		## samblaster mark duplicates for read id grouped reads. I do not coordinate sort the bam
		bowtie2 --threads 5  -X2000 -x {config[idx_bt2]} -1 {input[0]} -2 {input[1]} 2> {log.bowtie2} \
		| samblaster 2> {log.markdup} \
		| samtools view -Sb - > {output[0]}
		"""


# check number of reads mapped by samtools flagstat
rule flagstat_hg_bam:
    input:  "04aln/{case}_hg.sorted.bam"
    output: "04aln/{case}_hg.sorted.bam.flagstat"
    log:    "00log/{case}.flagstat_bam"
    threads: 1
    params: jobname = "{case}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

rule flagstat_toxo_bam:
    input:  "04aln/{case}_toxo.sorted.bam"
    output: "04aln/{case}_toxo.sorted.bam.flagstat"
    log:    "00log/{case}.flagstat_bam"
    threads: 1
    params: jobname = "{case}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

rule flagstat_control_bam:
    input:  "04aln/{control}.sorted.bam"
    output: "04aln/{control}.sorted.bam.flagstat"
    log:    "00log/{control}.flagstat_bam"
    threads: 1
    params: jobname = "{control}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """


rule ataqv:
	input: 
		ctl="04aln/{control}.sorted.bam",
		hum="04aln/{case}_hg.sorted.bam",
		
	output: "04aln/{control}.sorted.bam.ataqv.json","04aln/{case}_hg.sorted.bam.ataqv.json",
	log: "00log/{control}_ataqv.log","00log/{case}_ataqv.log"
	threads: 1
	params: jobname = "{input}"
	message: "ataqv quality control for {input}"
	shell:
		"""
		~/ataqv-1.0.0/bin/ataqv human {input.ctl} --metrics-file {output[0]} 2> {log}
		~/ataqv-1.0.0/bin/ataqv human {input.hum} --metrics-file {output[1]} 2> {log}
		"""

rule json_to_html:
	input: ATAQV_ALL
	output: "11ATAC_qc_html"
	log: "00log/ATAC_qc_html.log"
	threads: 1
	message: "compiling json files to html ATAC-seq QC"
	shell:
		"""
		#source activate root
		~/ataqv-1.0.0/bin/mkarv 11ATAC_qc_html {input}
		"""

    ## shifting the reads are only critical for TF footprint, for peak calling and making bigwigs, it should be fine using the bams without shifting
# https://sites.google.com/site/atacseqpublic/atac-seq-analysis-methods/offsetmethods
rule remove_chrM_bam:
	input: "04aln/{control}.sorted.bam"
	output: "06aln_exclude_chrM/{control}_exclude_chrM.sorted.bam", "06aln_exclude_chrM/{control}_exclude_chrM.sorted.bam.bai"
	log: "00log/{control}_exclude_chrM_bam.log"
	threads: 5
	message: "excluding chrM from bam {input} : {threads} threads"
	params: jobname = "{control}"
	shell:
		"""
		# remove duplicates and reads on chrM, coordinate sort the bam
		# samblaster expects name sorted bamq
		module add apps/samtools
		samtools view -h {input} | samblaster --removeDups \
		| grep -v -P '\tchrM\t' \
		| samtools view -Sb -F 4 - \
		| samtools sort -m 2G -@ 5 -T {input}.tmp -o {output[0]}
		samtools index {output[0]}
		"""

rule remove_chrM_bam_hg:
	input: "04aln/{case}_hg.sorted.bam"
	output: "06aln_exclude_chrM/{case}_hg_exclude_chrM.sorted.bam", "06aln_exclude_chrM/{case}_hg_exclude_chrM.sorted.bam.bai"
	log: "00log/{case}_exclude_chrM_bam.log"
	threads: 5
	message: "excluding chrM from bam {input} : {threads} threads"
	params: jobname = "{case}"
	shell:
		"""
		# remove duplicates and reads on chrM, coordinate sort the bam
		# samblaster expects name sorted bamq
		module add apps/samtools
		samtools view -h {input} | samblaster --removeDups \
		| grep -v -P '\tchrM\t' \
		| samtools view -Sb -F 4 - \
		| samtools sort -m 2G -@ 5 -T {input}.tmp -o {output[0]}
		samtools index {output[0]}
		"""
