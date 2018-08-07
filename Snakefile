import json
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output
import csv

# Globals ---------------------------------------------------------------------
configfile: "config.yml"

GENOME='~/toxo_kim/genomes/hg38.fa'
TOXO='~/toxo_kim/genomes/ToxoDB-38_TgondiiME49_Genome.fasta'
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


# Functions -------------------------------------------------------------------

def rstrip(text, suffix):
    # Remove a suffix from a string.
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

# Rules -----------------------------------------------------------------------

rule all:
    input:
        hum_cl=join(dirname(CASES),'clean_','hg38','.fq'),
        tox_cl=join(dirname(CASES),'clean_','ToxoDB-38_TgondiiME49_Genome','.fq')
rule clean_fastq:
    input:
        genome=GENOME,
        toxo=TOXO,
        fwd="/work/t/tekeller/atac_toxo/{cases}.1_val_1.fq.gz",
        rev="/work/t/tekeller/atac_toxo/{cases}.2_val_2.fq.gz"
    output:
       hum_cl=join(dirname(SAMPLES),'clean_','hg38','.fq'),
       tox_cl=join(dirname(SAMPLES),'clean_','ToxoDB-38_TgondiiME49_Genome','.fq')
    shell:
        """
        bbsplit.sh in={input.fwd} in2={input.rev} ref={input.genome},{input.toxo} basename={wildcards.sample}_clean_%.fq
        """

