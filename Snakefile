# **************************************
# * Snakefile for metassemble pipeline *
# **************************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input: "output/metaquast/report.html"

rule metaspades:
    input:
        r1 = config["path"]+"{sample}_1.fastq.gz",
        r2 = config["path"]+"{sample}_2.fastq.gz"
    output: "output/metaspades/{sample}/scaffolds.fasta"
    params:
        outdir = "output/metaspades/{sample}/",
        ec = "" if config["error_corr"] else "--only-assembler"
    conda: "metassemble_files/envs/metaspades_env.yaml"
    shell: "metaspades.py -o {params.outdir} -1 {input.r1} -2 {input.r2} -t 40 {params.ec}"

rule metaquast:
    input: expand("output/metaspades/{sample}/scaffolds.fasta", sample=SAMPLES)
    output: "output/metaquast/report.html"
    conda: "metassemble_files/envs/metaquast_env.yaml"
    shell: "metaquast.py {input} -o output/metaquast -t 40"

rule reformat:
    input:
        assembly = expand("output/metaspades/{sample}/scaffolds.fasta", sample=SAMPLES),
        r1 = config["path"]+"{sample}_1.fastq.gz",
        r2 = config["path"]+"{sample}_2.fastq.gz"
    output: expand("output/metaspades/{sample}/{sample}_idxstats.txt", sample=SAMPLES)
    params:
        s = "{sample}"
    conda: "metassemble_files/envs/metassemble_env.yaml"
    shell:
            "bowtie2-build --seed 1 {input.assembly} output/metaspades/{params.s}/bowtie_db/{params.s}_db; bowtie2 --sensitive-local -p 40 --seed 1 -x output/metaspades/{params.s}/bowtie_db/{params.s}_db -1 {input.r1} -2 {input.r2} -S output/metaspades/{params.s}/{params.s}.sam; samtools faidx output/metaspades/{params.s}/scaffolds.fasta; samtools import output/metaspades/{params.s}/scaffolds.fai output/metaspades/{params.s}/{params.s}.sam output/metaspades/{params.s}/{params.s}.bam; samtools sort output/metaspades/{params.s}/{params.s}.bam -o output/metaspades/{params.s}/{params.s}_sorted.bam; samtools index output/metaspades/{params.s}/{params.s}_sorted.bam; samtools idxstats output/metaspades/{params.s}/{params.s}_sorted.bam > output/metaspades/{params.s}/{params.s}_idxstats.txt"

rule get_count_table:
    input: expand("output/metaspades/{sample}/{sample}_idxstats.txt", sample=SAMPLES)
    output: expand("output/metaspades/{sample}/{sample}_counts.txt", sample=SAMPLES)
    conda: "metassemble_files/envs/python27_env.yaml"
    shell:
            "python metassemble_files/scripts/get_count_table.py {input} > {output}"
