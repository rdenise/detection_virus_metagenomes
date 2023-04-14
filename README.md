# Snakemake workflow: Virome pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.24.2-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/617537743.svg)](https://zenodo.org/badge/latestdoi/617537743)

## Aim

The workflow takes raw reads from metagenomes as input, performs quality control and trimming using FastQC and FastP, and then assembles the reads using MEGAHIT and SPAdes. The resulting contigs are then analyzed using viral identification tool (geNomad) to identify putative viral sequences.

Next, the identified viral sequences are taxonomically classified using genomad and clusterisation using ANI, and potential hosts for the viruses are inferred using iPhoP. Finally, infer their lifestyle using Baphlip. 

The output of each step is stored in a separate directory, and the workflow is managed using Snakemake to ensure efficient use of computing resources and reproducibility of results 

## Installation

### Step 1: install Snakemake and Snakedeploy

Snakemake and Snakedeploy are best installed via the [Mamba package manager](https://github.com/mamba-org/mamba) (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge). For other options see [here](https://github.com/mamba-org/mamba).

To install Mamba by conda, run

```shell
conda install mamba -n base -c conda-forge
```

Given that Mamba is installed, run 

```shell
mamba create -c bioconda -c conda-forge --name snakemake snakemake snakedeploy
```

to install both Snakemake and Snakedeploy in an isolated environment. 

#### Notes 

For all following commands (step 2 to 5) ensure that your snakemake environment is activated via 

```shell
conda activate snakemake
```

### Step 2: deploy workflow

Given that Snakemake and Snakedeploy are installed and available (see Step 1), the workflow can be deployed as follows.

First, create an appropriate project working directory on your system in the place of your choice as follow (note to change the path and file name to the one you want to create): : 

```shell
mkdir path/to/project-workdir
```

Then go to your project working directory as follow:

```shell
cd path/to/project-workdir
```

In all following steps, we will assume that you are inside of that directory.

Second, run 

```shell
snakedeploy deploy-workflow https://github.com/rdenise/detection_virus_metagenomes . --tag 0.0.2
```

Snakedeploy will create two folders `workflow` and `config`. The former contains the deployment of the chosen workflow as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows), the latter contains configuration files which will be modified in the next step in order to configure the workflow to your needs. Later, when executing the workflow, Snakemake will automatically find the main `Snakefile` in the `workflow` subfolder.

### Step 3: configure workflow

#### General settings

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.  

Missing values can be specified by empty columns or by writing `NA`.

### Step 4: run the workflow

Given that the workflow has been properly deployed and configured, it can be executed as follows.

For running the workflow while deploying any necessary software via conda, run Snakemake with 

```shell
snakemake --cores 10 --use-conda 
```

### Step 5: MutiQC report

After finalizing your data analysis, a report is produce by MultiQC summariing some of the key statistics

## Walk-Through and File Production

This workflow consists of steps called rules that take input files and create output files. Here is a description of the pipeline.

1. As snakemake is set up, there is a last rule, called `all`, that serves to call the last output files and make sure they were created.

2. A folder containing your work will be created:

```
[output_name]                          <- Main results folder
├── databases                          <- Folder containing databases used in the analysis
│   ├── checkv_db                      <- CheckV database for viral genome completeness and contamination assessment
│   ├── contigs                        <- Folder containing contig FASTA files
│   └── reads_trimmed                  <- Folder containing trimmed read FASTQ files
├── logs                               <- Folder containing log files for each analysis step
├── processed_files                    <- Folder containing processed files resulting from the analysis
│   ├── assemblies                     <- Folder containing assembled contigs
│   ├── blast                          <- Folder containing BLAST output files
│   ├── bowtie2                        <- Folder containing Bowtie2 output files
│   ├── checkv                         <- Folder containing CheckV output files
│   ├── genomad                        <- Folder containing downloaded GenomAD data
│   ├── otu                            <- Folder containing OTU clustering output files
│   └── samtools                       <- Folder containing Samtools output files
├── qc                                 <- Folder containing quality control reports
│   ├── fastqc                         <- FastQC report files for each input read file
│   ├── multiqc_report.trimmed_data    <- MultiQC report for trimmed reads
│   └── multiqc_report.untrimmed_data  <- MultiQC report for untrimmed reads
└── results                            <- Folder containing final analysis results
    ├── bacphlip_out                   <- Output files for BacPhlip viral protein prediction tool
    ├── iphop                          <- Output files for IPHOP prophage prediction tool
    ├── taxonomy                       <- Folder containing taxonomy assignment output files
    └── viral_contigs                  <- Folder containing viral contigs identified in the analysis
```

3. When restarting the pipeline, the software will check if you made any changes in the seed file before running. If changes have been made, it will run what is necessary, else nothing will happen.

### Pipeline in image 

#### Normal behavior

<p align="center">
  <img src="doc/dag.png?raw=true" height="400">
</p>

