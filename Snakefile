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
    input:
        "output/metaquast/report.html",
        "output/counts.txt"

rule metaspades:
    input:
        r1 = config["forward"],
        r2 = config["reverse"]
    output: "output/assembly/scaffolds.fasta"
    params:
        ec = "" if config["error_corr"] else "--only-assembler"
    conda: "metassemble_files/envs/metaspades_env.yaml"
    shell: "metaspades.py -o output/assembly/ -1 {input.r1} -2 {input.r2} -t 56 {params.ec}"

rule metaquast:
    input: "output/assembly/scaffolds.fasta"
    output: "output/metaquast/report.html"
    conda: "metassemble_files/envs/metaquast_env.yaml"
    shell: "metaquast.py {input} -o output/metaquast -t 28"

rule index_reference:
    input: "output/assembly/scaffolds.fasta"
    output: "output/assembly/scaffolds.fasta.fai"
    conda: "metassemble_files/envs/metassemble_env.yaml"
    shell: "bowtie2-build --seed 1 {input} output/assembly/bowtie_db/reference; samtools faidx {input}"

rule map:
    input:
        r1 = config["path"]+"{sample}_bmtagged_1.fastq",
        r2 = config["path"]+"{sample}_bmtagged_2.fastq",
        ref = "output/assembly/scaffolds.fasta.fai"
    output: "output/mapping/{sample}/{sample}_idxstats.txt"
    params:
        s = "{sample}"
    conda: "metassemble_files/envs/metassemble_env.yaml"
    shell:
            "bowtie2 --sensitive-local -p 40 --seed 1 -x output/assembly/bowtie_db/reference -1 {input.r1} -2 {input.r2} -S output/mapping/{params.s}/{params.s}.sam; samtools import output/assembly/scaffolds.fasta.fai output/mapping/{params.s}/{params.s}.sam output/mapping/{params.s}/{params.s}.bam; samtools sort output/mapping/{params.s}/{params.s}.bam -o output/mapping/{params.s}/{params.s}_sorted.bam; samtools index output/mapping/{params.s}/{params.s}_sorted.bam; samtools idxstats output/mapping/{params.s}/{params.s}_sorted.bam > {output}"

rule summarize_counts:
    input: expand("output/mapping/{sample}/{sample}_idxstats.txt", sample=SAMPLES)
    output: "output/counts.txt"
    conda: "metassemble_files/envs/python27_env.yaml"
    shell:
            "python metassemble_files/scripts/get_count_table.py output/mapping/*/*_idxstats.txt > {output}"
