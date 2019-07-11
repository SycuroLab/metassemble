# **********************************
# * Snakefile for metassemble pipeline *
# **********************************

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
        outdir = "output/metaspades/{sample}/"
        ec = "" if config["error_corr"] else "--only-assembler"
    conda: "metassemble_files/envs/metaspades_env.yaml"
    shell: "metaspades.py -o {params.outdir} -1 {input.r1} -2 {input.r2} -t 40 {params.ec}"

rule metaquast:
    input: expand("output/metaspades/{sample}/scaffolds.fasta", sample=SAMPLES)
    output: "output/metaquast/report.html"
    conda: "metassemble_files/envs/metaquast_env.yaml"
    shell: "metaquast.py {input} -o output/metaquast -t 40"
