# **********************************
# * Snakefile for metassemble pipeline *
# **********************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

import subprocess
from os.path import join

# **** Rules ****

rule all:
    input:
        # List final outputs here

rule metaspades:
    input:
        r1 = join(config["path"], "{sample}_1.fastq.gz"),
        r2 = join(config["path"], "{sample}_2.fastq.gz")
    output:

    conda: "metassemble_files/envs/metaspades_env.yaml"
    shell: "spades.py --meta -o data/metaspades/"

rule megahit:
    input:
        r1 = join(config["path"], "{sample}_1.fastq.gz"),
        r2 = join(config["path"], "{sample}_2.fastq.gz")
    output
    conda: "metassemble_files/envs/megahit_env.yaml"
    shell: "megahit usage"

rule metaquast:
    input:

    output:

    conda: "metassemble_files/envs/metaquast_env.yaml"
    shell:
            "super duper long usage of metaquast here "
            "second line of super long usage command "
            "maybe even this long whoa"

rule checkm:
    input:

    output:

    conda:
    shell:
