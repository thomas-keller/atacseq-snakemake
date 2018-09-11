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


CASE_FILES=expand("/work/t/tekeller/atac_toxo/{case}.1_val_1.fq.gz")

logdir = os.path.join(os.getcwd(), "logs/slurm")
os.makedirs(logdir, exist_ok=True)

# Functions -------------------------------------------------------------------

def rstrip(text, suffix):
    # Remove a suffix from a string.
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

# Rules -----------------------------------------------------------------------

CONTROL_MERGED_FASTQ = expand("02seq/{control}_merged.fq.gz", control = CONTROLS)
CASE_CLEAN_HG=expand("/work/t/tekeller/atac_toxo/01cln/{case}_clean_hg38.fq",case=CASES)
CASE_CLEAN_TOXO=expand("/work/t/tekeller/atac_toxo/01cln/{case}_clean_ToxoDB-38_TgondiiME49_Genome.fq",case=CASES)
ALL_SAMPLES=expand("work/t/tekeller/atac_toxo/{sample}.1_val_1.fq.gz",sample=SAMPLES)

 ATAQV_CTL=expand("04aln/{control}.sorted.bam.ataqv.json",control=CONTROLS)
 ATAQV_HG=expand("04aln/{case}_hg.sorted.bam.ataqv.json",case=CASES)
ATAQV_ALL=ATAQV_CTL+ATAQV_HG

ALN_HG=expand("04aln/{case}_hg.sorted.bam",case=CASES)
ALN_TOXO=expand("04aln/{case}_toxo.sorted.bam",case=CASES)
ALN_CTL=expand("04aln/{control}_toxo.sorted.bam",control=CONTROLS)
ALN_ALL=ALN_HG+ALN_TOXO+ALN_CTL

FLAG_HG=expand("04aln/{case}_hg.sorted.bam.flagstat",case=CASES)
FLAG_TOXO=expand("04aln/{case}_toxo.sorted.bam.flagstat",case=CASES)
FLAG_CTL=expand("04aln/{control}_toxo.sorted.bam.flagstat",control=CONTROLS)
FLAG_ALL=FLAG_HG+FLAG_TOXO+FLAG_CTL

NUCL_CASE=expand("09nucleoATAC/{case}_nucleoATAC.occpeaks.bed.gz",case=CASES)
NUCL_CTL=expand("09nucleoATAC/{control}_nucleoATAC.occpeaks.bed.gz",control=CONTROLS)
NUCL_ALL=NUCL_CASE+NUCL_CTL

print(CASES)
rule all:
    input:
        CONTROL_MERGED_FASTQ + CASE_CLEAN_HG + CASE_CLEAN_TOXO+ALN_ALL+FLAG_ALL+NUCL_ALL

CONTROL_FILES=expand('/work/t/tekeller/atac_toxo/{control}.1_val_1.fq.gz',control=CONTROLS)


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
		mv /work/t/tekeller/atac_toxo/{wildcards.sample}.1_val_1.trimmed.fq.gz 01trim/
		mv /work/t/tekeller/atac_toxo/{wildcards.sample}.2_val_2.trimmed.fq.gz 01trim/
 		"""


rule clean_fastq:
    input:
        genome=GENOME,
        toxo=TOXO,
        fwd="02trim/{case}.1_val_1.fq.gz",
        rev="02trim/{case}.2_val_2.fq.gz"
    output:
       hum_cl="/work/t/tekeller/atac_toxo/03cln/{case}_clean_hg38.fq",
       tox_cl="/work/t/tekeller/atac_toxo/03cln/{case}_clean_ToxoDB-38_TgondiiME49_Genome.fq"
    shell:
        """
        bbsplit.sh -Xmx28g -t=8 in={input.fwd} in2={input.rev} basename=02cln/{wildcards.case}_clean_%.fq
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
		bowtie2 --threads 5  -X2000 -x {config[idx_bt2]} -1 {input[0]} -2 {input[1]} 2> {log.bowtie2} \
		| samblaster 2> {log.markdup} \
		| samtools view -Sb - > {output[0]}
		"""

rule align_cases_toxo:
	input: "/work/t/tekeller/atac_toxo/03cln/{case}_clean_ToxoDB-38_TgondiiME49_Genome.fq"
	output: "04aln/{case}_toxo.sorted.bam", "00log/{case}.align"
	params: jobname = "{sample}"
	message: "aligning {input}: {threads} threads"
	log:
		bowtie2 = "00log/{case}_toxo.align",
		markdup = "00log/{case}_toxo.markdup"
	shell:
		"""
		## samblaster mark duplicates for read id grouped reads. I do not coordinate sort the bam
		bowtie2 --threads 5  -X2000 -x {config[idx_bt2]} -1 {input[0]} -2 {input[1]} 2> {log.bowtie2} \
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
	log: "00log/{sample}_ataqv.log"
	threads: 1
	params: jobname = "{sample}"
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
	output: "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam", "06aln_exclude_chrM/{sample}_exclude_chrM.sorted.bam.bai"
	log: "00log/{sample}_exclude_chrM_bam.log"
	threads: 5
	message: "excluding chrM from bam {input} : {threads} threads"
	params: jobname = "{sample}"
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
	params: jobname = "{sample}"
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

rule phantom_peak_qual:
    input: "06aln_exclude_chrM/{control}_exclude_chrM.sorted.bam"
    output: "05phantompeakqual/{control}_phantom.txt"
    log: "00log/{control}_phantompeakqual.log"
    threads: 4
    params: jobname = "{control}"
    message: "phantompeakqual for {input} : {threads} threads"
    shell:
        """
		source activate bioconductor
		module add apps/samtools
        Rscript  ~/miniconda3/envs/bioconductor/bin/run_spp_nodups.R -c={input} -savp -rf  -p=4 -odir=05phantompeakqual  -out={output} -tmpdir=05phantompeakqual 2> {log}
        """

rule phantom_peak_qual_hg:
    input: "06aln_exclude_chrM/{case}_hg_exclude_chrM.sorted.bam"
    output: "05phantompeakqual/{case}_hg_phantom.txt"
    log: "00log/{case}_phantompeakqual.log"
    threads: 4
    params: jobname = "{case}"
    message: "phantompeakqual for {input} : {threads} threads"
    shell:
        """
		source activate bioconductor
		module add apps/samtools
        Rscript  ~/miniconda3/envs/bioconductor/bin/run_spp_nodups.R -c={input} -savp -rf  -p=4 -odir=05phantompeakqual  -out={output} -tmpdir=05phantompeakqual 2> {log}
        """

## consider how to reuse the rules.
rule flagstat_chrM_exclude_bam_case:
    input:  "06aln_exclude_chrM/{control}_exclude_chrM.sorted.bam"
    output: "06aln_exclude_chrM/{control}_exclude_chrM.sorted.bam.flagstat"
    log:    "00log/{control}_exclude_chrM_flagstat_bam"
    threads: 1
    params: jobname = "{control}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

rule flagstat_chrM_exclude_bam_hg:
    input:  "06aln_exclude_chrM/{case}_exclude_chrM.sorted.bam"
    output: "06aln_exclude_chrM/{case}_exclude_chrM.sorted.bam.flagstat"
    log:    "00log/{sample}_exclude_chrM_flagstat_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

rule down_sample:
    input: "06aln_exclude_chrM/{control}_exclude_chrM.sorted.bam", "06aln_exclude_chrM/{control}_exclude_chrM.sorted.bam.bai",
	 	   "06aln_exclude_chrM/{control}_exclude_chrM.sorted.bam.flagstat"
    output: "06aln_downsample/{control}-downsample.sorted.bam", "06aln_downsample/{control}-downsample.sorted.bam.bai"
    log: "00log/{sample}_downsample.log"
    threads: 5
    params: jobname = "{control}"
    message: "downsampling for {input}"
    run:
        import re
        import subprocess
        with open (input[2], "r") as f:
            # fifth line contains the number of mapped reads
            line = f.readlines()[5]
            match_number = re.match(r'(\d.+) \+.+', line)
			## how many paired reads, roughly total #reads/2
            total_reads = float(match_number.group(1))/2

        target_reads = config["target_reads"] # 15million reads  by default, set up in the config.yaml file
        if total_reads > target_reads:
            down_rate = target_reads/total_reads
        else:
            down_rate = 1

        shell("sambamba view -f bam -t 5 --subsampling-seed=3 -s {rate} {inbam} | samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} 2> {log}".format(rate = down_rate, inbam = input[0], outbam = output[0], log = log))

        shell("samtools index {outbam}".format(outbam = output[0]))

rule down_sample_hg:
    input: "06aln_exclude_chrM/{case}_hg_exclude_chrM.sorted.bam", "06aln_exclude_chrM/{case}_hg_exclude_chrM.sorted.bam.bai",
	 	   "06aln_exclude_chrM/{case}_exclude_chrM.sorted.bam.flagstat"
    output: "06aln_downsample/{case}-hg-downsample.sorted.bam", "06aln_downsample/{case}-hg-downsample.sorted.bam.bai"
    log: "00log/{case}_downsample.log"
    threads: 5
    params: jobname = "{case}"
    message: "downsampling for {input}"
    run:
        import re
        import subprocess
        with open (input[2], "r") as f:
            # fifth line contains the number of mapped reads
            line = f.readlines()[5]
            match_number = re.match(r'(\d.+) \+.+', line)
			## how many paired reads, roughly total #reads/2
            total_reads = float(match_number.group(1))/2

        target_reads = config["target_reads"] # 15million reads  by default, set up in the config.yaml file
        if total_reads > target_reads:
            down_rate = target_reads/total_reads
        else:
            down_rate = 1

        shell("sambamba view -f bam -t 5 --subsampling-seed=3 -s {rate} {inbam} | samtools sort -m 2G -@ 5 -T {outbam}.tmp > {outbam} 2> {log}".format(rate = down_rate, inbam = input[0], outbam = output[0], log = log))

        shell("samtools index {outbam}".format(outbam = output[0]))


rule make_bigwigs:
    input : "06aln_downsample/{control}-downsample.sorted.bam", "06aln_downsample/{control}-downsample.sorted.bam.bai"
    output: "07bigwig/{control}.bw"
    log: "00log/{control}.makebw"
    threads: 5
    params: jobname = "{control}"
    message: "making bigwig for {input} : {threads} threads"
    shell:
        """
    	#source activate root
    	# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
        bamCoverage -b {input[0]} --ignoreDuplicates --skipNonCoveredRegions --normalizeUsingRPKM -p 5 --extendReads -o {output} 2> {log}
        """

rule make_bigwigs_hg:
    input : "06aln_downsample/{case}-hg-downsample.sorted.bam", "06aln_downsample/{case}-hg-downsample.sorted.bam.bai"
    output: "07bigwig/{case}_hg.bw"
    log: "00log/{case}.makebw"
    threads: 5
    params: jobname = "{case}"
    message: "making bigwig for {input} : {threads} threads"
    shell:
        """
    	#source activate snakemake
    	# no window smoothing is done, for paired-end, bamCoverage will extend the length to the fragement length of the paired reads
        bamCoverage -b {input[0]} --ignoreDuplicates --skipNonCoveredRegions --normalizeUsingRPKM -p 5 --extendReads -o {output} 2> {log}
        """
# https://github.com/taoliu/MACS/issues/145
rule call_peaks_macs2:
    input: "06aln_downsample/{control}-downsample.sorted.bam", "06aln_downsample/{control}-downsample.sorted.bam.bai"
    output: bed = "08peak_macs2/{control}_macs2_peaks.broadPeak"
    log: "00log/{control}_call_peaks_macs2.log"
    params:
        name = "{control}_macs2",
        jobname = "{control}"
    message: "call_peaks macs2 {input}: {threads} threads"
    shell:
        """
       source activate macs2
       ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input[0]} \
            --keep-dup all -f BAMPE -g {config[macs2_g]} \
            --outdir 08peak_macs2 -n {params.name} -p {config[macs2_pvalue]} \
            --broad --broad-cutoff {config[macs2_pvalue_broad]} &> {log}
        """

rule call_peaks_macs2_hg:
    input: "06aln_downsample/{case}-hg-downsample.sorted.bam", "06aln_downsample/{case}-hg-downsample.sorted.bam.bai"
    output: bed = "08peak_macs2/{case}_hg_macs2_peaks.broadPeak"
    log: "00log/{case}_call_peaks_macs2.log"
    params:
        name = "{case}_macs2",
        jobname = "{case}"
    message: "call_peaks macs2 {input}: {threads} threads"
    shell:
        """
       source activate macs2
       ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input[0]} \
            --keep-dup all -f BAMPE -g {config[macs2_g]} \
            --outdir 08peak_macs2 -n {params.name} -p {config[macs2_pvalue]} \
            --broad --broad-cutoff {config[macs2_pvalue_broad]} &> {log}
        """



rule multiQC:
    input :
        expand("00log/{case}_hg.align", case = CASES),
		expand("00log/{case}_toxo.align", case = CASES),
		expand("00log/{control}.align", case = CASES1),
        expand("04aln/{case}_hg.sorted.bam.flagstat", case= CASES),
		 expand("04aln/{control}_toxo.sorted.bam.flagstat", case= CASES),
        expand("02fqc/{case}.{read}.fastqc.zip", sample = ALL_SAMPLES, read = ["1_val_1", "2_val_2"])
    output: "10multiQC/multiQC_log.html"
    log: "00log/multiqc.log"
    message: "multiqc for all logs"
    shell:
        """
        multiqc 02fqc 04aln 00log -o 10multiQC -d -f -v -n multiQC_log 2> {log}
        """



## extend the broad peak a bit for nucelosome analysis by nuceloATAC
rule make_bed_nucleoATAC_ctl:
	input: "08peak_macs2/{control}_hg_macs2_peaks.broadPeak"
	output: "09nucleoATAC/{control}_hg_nucleo.bed"
	log: "00log/{control}_make_nucleoATAC_bed.log"
	threads: 1
	message: "making bed for nucleoATAC from {input}"
	params: jobname= "{control}"
	shell:
		"""
		cat {input} | bedtools slop -b 200 -g {config[genome_size]} | sort -k1,1 -k2,2n | bedtools merge > {output} 2> {log}
		"""

rule make_bed_nucleoATAC_hg:
	input: "08peak_macs2/{case}_hg_macs2_peaks.broadPeak"
	output: "09nucleoATAC/{case}_hg_nucleo.bed"
	log: "00log/{case}_make_nucleoATAC_bed.log"
	threads: 1
	message: "making bed for nucleoATAC from {input}"
	params: jobname= "{case}"
	shell:
		"""
		cat {input} | bedtools slop -b 200 -g {config[genome_size]} | sort -k1,1 -k2,2n | bedtools merge > {output} 2> {log}
		"""

## nucleoATAC works on the non-shifted bam, and shift the reads internally!
# https://github.com/GreenleafLab/NucleoATAC/issues/58
rule nucleoATAC_ctl:
	input: "06aln_downsample/{control}-downsample.sorted.bam", "06aln_downsample/{control}-downsample.sorted.bam.bai", "09nucleoATAC/{control}_hg_nucleo.bed"
	output: "09nucleoATAC/{control}_nucleoATAC.occpeaks.bed.gz"
	log: "00log/{control}_nucleoATAC.log"
	threads: 5
	message: "calling nucleosome by nucleoATAC for {input} : {threads} threads"
	params:
		jobname = "{control}",
		outputdir = os.path.dirname(srcdir("00log"))
	shell:
		"""
		#source activate root
		cd 09nucleoATAC
		nucleoatac run --bed {params.outputdir}/{input[2]} --bam {params.outputdir}/{input[0]} --cores 5 --fasta {config[genome_fasta]} --out {wildcards.sample} 2> {params.outputdir}/{log}
		"""

rule nucleoATAC_hg:
	input: "06aln_downsample/{case}-hg-downsample.sorted.bam", "06aln_downsample/{case}-hg-downsample.sorted.bam.bai", "09nucleoATAC/{case}_hg_nucleo.bed"
	output: "09nucleoATAC/{case}_nucleoATAC.occpeaks.bed.gz"
	log: "00log/{case}_nucleoATAC.log"
	threads: 5
	message: "calling nucleosome by nucleoATAC for {input} : {threads} threads"
	params:
		jobname = "{case}",
		outputdir = os.path.dirname(srcdir("00log"))
	shell:
		"""
		#source activate root
		cd 09nucleoATAC
		nucleoatac run --bed {params.outputdir}/{input[2]} --bam {params.outputdir}/{input[0]} --cores 5 --fasta {config[genome_fasta]} --out {wildcards.sample} 2> {params.outputdir}/{log}
		"""