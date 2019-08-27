# metassemble

Bioinformatics pipeline for assembly of shotgun metagenomic data, using [metaSPAdes](https://www.ncbi.nlm.nih.gov/pubmed/28298430) and [MetaQUAST](https://www.ncbi.nlm.nih.gov/pubmed/26614127).

## Overview

This pipeline is written in snakemake and designed to automate and control the submission of processes to the Synergy server at the University of Calgary. Developed by Alana Schick for the lab of Dr. Laura Sycuro. 

Input: paired-end fastq files.

Output: assembled fasta files and a detailed report of the quality of the assemblies.

## Installation

To use this pipeline, navigate to your project directory and clone this repository into that directory using the following command:

```
git clone https://github.com/alanaschick/metassemble.git metassemble
```

Note: you need to have **conda** and **snakemake** installed in order to run this. To install conda, see the instructions [here](https://github.com/ucvm/synergy/wiki). 

To install snakemake using conda, run the following line:

```
conda install -c bioconda -c conda-forge snakemake
```

See the snakemake installation [webpage](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for further details.

## Config file

All the parameters required to run this pipeline are specified in a config file, written in yaml. See/modify the provided example file with your custom parameters, called `config.yaml`. This is the only file that should be modified before running the pipeline. Make sure to follow the syntax in the example file in terms of when to use quotations around parameters.

## Data and list of files

Specify the full path to the directory that contains your data files in the config file. You also need to have a list of sample names which contains the names of the samples to run the pipeline on, one sample per line. You can run this pipeline on any number or subset of your samples. Sample names should include everything up to the R1/R2 (or 1/2) part of the file names of the raw fastq files. Specify the path and name of your list in the config file.

If there are many samples, it may be convenient to generate the list of files using the following command, replacing `R1_001.fastq.gz` with the general suffix of your files:

```
ls | grep R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' > list_files.txt
```

## Description of other parameters
| Parameter | Description |
| -------------- | --------------- |
| list_files | Full path and name of your sample list. |
| path | Location of input files. |
| error_corr | If TRUE (default), include read error correction during metaspades run. |

## Running the pipeline on Synergy

Test the pipeline by running `snakemake -np`. This command prints out the commands to be run without actually running them. 

To run the pipeline on the Synergy compute cluster, enter the following command from the project directory:

```
snakemake --cluster-config cluster.json --cluster 'bsub -n {cluster.n} -R {cluster.resources} -W {cluster.walllim} -We {cluster.time} -M {cluster.maxmem} -oo {cluster.output} -e {cluster.error}' --jobs 100 --use-conda
```
The above command submits jobs to Synergy, one for each sample and step of the QC pipeline. Note: the file `cluster.json` contains the parameters for the LSF job submission system that Synergy uses. In most cases, this file should not be modified.

## Results and log files

Snakemake will create a directory for the results of the pipeline as well as a directory for log files. Log files of each step of the pipeline will be written to the `logs` directory.

## Pipeline summary

![Rulegraph](./metassemble_files/metassemble_dag.png)

### Steps

1) Metagenome assembly using metaSPAdes. Includes error correction by default. To disable, this, set the `error_corr` parameter to `FALSE` in the `config.yaml` file. See [paper](https://genome.cshlp.org/content/27/5/824.short) for details about metaSPAdes.

2) QC on metagenome assembly MetaQUAST. More details [here](http://quast.sourceforge.net/metaquast).


