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
CASE_CLEAN_HG=expand("/work/t/tekeller/atac_toxo/01cln/{case}_clean_hg38.fq",case=CASES)
CASE_CLEAN_TOXO=expand("/work/t/tekeller/atac_toxo/01cln/{case}_clean_ToxoDB-38_TgondiiME49_Genome.fq",case=CASES)
ALL_SAMPLES=expand("work/t/tekeller/atac_toxo/{sample}.{read}.fq.gz",sample=SAMPLES,read=['1_val_1','2_val_2'])

FASTQC_ALL=expand('/work/t/tekeller/atac_toxo/01fqc/{sample}.{read}.fastqc.zip',sample=SAMPLES,read=['1_val_1','2_val_2'])
# ATAQV_CTL=expand("04aln/{control}.sorted.bam.ataqv.json",control=CONTROLS)
# ATAQV_HG=expand("04aln/{case}_hg.sorted.bam.ataqv.json",case=CASES)
#ATAQV_ALL=ATAQV_CTL+ATAQV_HG

#ALN_HG=expand("04aln/{case}_hg.sorted.bam",case=CASES)
#ALN_TOXO=expand("04aln/{case}_toxo.sorted.bam",case=CASES)
#ALN_ALL=ALN_HG+ALN_TOXO+ALN_CTL
#ALN_CTL=expand("04aln/{control}.sorted.bam",control=CONTROLS)

#FLAG_HG=expand("04aln/{case}_hg.sorted.bam.flagstat",case=CASES)
#FLAG_TOXO=expand("04aln/{case}_toxo.sorted.bam.flagstat",case=CASES)
#FLAG_CTL=expand("04aln/{control}_toxo.sorted.bam.flagstat",control=CONTROLS)
#FLAG_ALL=FLAG_HG+FLAG_TOXO+FLAG_CTL

#NUCL_CASE=expand("09nucleoATAC/{case}_nucleoATAC.occpeaks.bed.gz",case=CASES)
#NUCL_CTL=expand("09nucleoATAC/{control}_nucleoATAC.occpeaks.bed.gz",control=CONTROLS)
#NUCL_ALL=NUCL_CASE+NUCL_CTL

print(CASES)
print(CONTROLS)
rule all:
    input:
        #CONTROL_MERGED_FASTQ + CASE_CLEAN_HG + CASE_CLEAN_TOXO+ALN_ALL+FLAG_ALL+NUCL_ALL
        ALL_SAMPLES+ CASE_CLEAN_HG + CASE_CLEAN_TOXO

CONTROL_FILES=expand('/work/t/tekeller/atac_toxo/{control}.1_val_1.fq.gz',control=CONTROLS)
ALL_SAMPLES=expand("work/t/tekeller/atac_toxo/{sample}.1_val_1.fq.gz",sample=SAMPLES)
ALL_CASES=expand("work/t/tekeller/atac_toxo/{case}.1_val_1.fq.gz",case=CASES)

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
		mv /work/t/tekeller/atac_toxo/{wildcards.sample}.1_val_1.trimmed.fq.gz 02trim/
		mv /work/t/tekeller/atac_toxo/{wildcards.sample}.2_val_2.trimmed.fq.gz 02trim/
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
        bbsplit.sh -Xmx28g -t=8 in={input.fwd} in2={input.rev} out_hg38.fq={output.hum_cl} out_ToxoDB-38_TgondiiME49_Genome.fq={output.tox_cl}
        """

