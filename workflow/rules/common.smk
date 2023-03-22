##########################################################################
##########################################################################
##
##                                Library
##
##########################################################################
##########################################################################

import os, sys
import pandas as pd
import numpy as np
from snakemake.utils import validate
import glob

##########################################################################
##########################################################################
##
##                               Functions
##
##########################################################################
##########################################################################


def get_final_output(outdir, reads=False):
    """
    Generate final output name
    """
    final_output = []

    if reads:
        # multiqc untrimmed
        final_output += (
            os.path.join(
                outdir,
                "qc",
                "multiqc_report.untrimmed.html"
            ),
            os.path.join(
                outdir,
                "qc",
                "multiqc_report.trimmed.html",
            ),
        )

    # Host annotation
    final_output += (
        os.path.join(
            outdir,
            "results",
            "iphop",
        ),
    )

    # Lifestyle annotation
    final_output += (
        os.path.join(
            outdir,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna.bacphlip",
        )
    )

    # Taxonomy annotation
    final_output += (
        os.path.join(
            outdir,
            "results",
            "taxonomy",
            "blast_table_species.tsv",
        ),
    )

    # Taxonomy annotation
    final_output += (
            os.path.join(
            outdir,
            "processed_files",
            "genomad",
            "viral_contigs",
            "all_contigs.over3kb.nr_annotate",
            "all_contigs.over3kb.nr_taxonomy.tsv",
        ),
    )

    # Taxonomy annotation
    final_output += (
        os.path.join(
            outdir,
            "processed_files",
            "otu",
            "final_vOTUs_representative_contigs.txt",
        ),
    )

    return final_output


##########################################################################


def create_folder(mypath):
    """
    Created the folder that I need to store my result if it doesn"t exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return


##########################################################################


def max_len_seq(file_fasta, ext_compress):
    max_len = 0
    tmp_len = 0
    if ext_compress == "tar.gz":
        import tarfile

        with tarfile.open(file_fasta, "r:gz") as tar:
            for tarinfo in tar:
                f = tar.extractfile(tarinfo.name)
                # To get the str instead of bytes str
                # Decode with proper coding, e.g. utf-8
                content = f.read().decode("utf-8", errors="ignore")
                # Split the long str into lines
                # Specify your line-sep: e.g. \n
                lines = content.split("\n")

                for line in lines:
                    if line.startswith(">"):
                        max_len = max(max_len, tmp_len)
                        tmp_len = 0
                    else:
                        tmp_len += len(line)
    elif ext_compress == "gz":
        import gzip

        with gzip.open(file_fasta, "rt") as r_file:
            for line in r_file:
                if line.startswith(">"):
                    max_len = max(max_len, tmp_len)
                    tmp_len = 0
                else:
                    tmp_len += len(line)

    elif ext_compress == "":
        with open(file_fasta, "rt") as r_file:
            for line in r_file:
                if line.startswith(">"):
                    max_len = max(max_len, tmp_len)
                    tmp_len = 0
                else:
                    tmp_len += len(line)

    return max_len


##########################################################################
##########################################################################
##
##                                Variables
##
##########################################################################
##########################################################################

# Validation of the config.yaml file
validate(config, schema="../schemas/config.schema.yaml")

##########################################################################
##########################################################################
##
##                        Core configuration
##
##########################################################################
##########################################################################

## Store some workflow metadata
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_basedir_short__"] = os.path.basename(workflow.basedir)
config["__workflow_workdir__"] = os.getcwd()

if workflow.config_args:
    tmp_config_arg = '" '.join(workflow.config_args).replace("=", '="')
    config["__config_args__"] = f' -C {tmp_config_arg}"'
else:
    config["__config_args__"] = ""

with open(os.path.join(workflow.basedir, "../config/VERSION"), "rt") as version:
    url = "https://github.com/rdenise/viral_detection/releases/tag"
    config["__workflow_version__"] = version.readline()
    config["__workflow_version_link__"] = f"{url}/{config['__workflow_version__']}"


##########################################################################
##########################################################################
##
##                           Options
##
##########################################################################
##########################################################################

# Result folder
OUTPUT_FOLDER = config["output_folder"]
# Adding to config for report
config["__output_folder__"] = os.path.abspath(OUTPUT_FOLDER)

# Options for blastn
blast_evalue = config["default_blast_option"]["e_val"]
blast_coverage = config["default_blast_option"]["coverage"]
blast_pident = config["default_blast_option"]["pident"]

blast_database = config["default_blast_option"]["nt"]

DB_DICT = config["databases"]

# path to contigs sheet (TSV format, columns: contig_name, path_contig)
CONTIGS_FOLDER = config["metagenomes"]["assemble_contigs"]

if not config["metagenomes"]["contigs_ext"].startswith("."):
    CONTIGS_EXT = f".{config['metagenomes']['contigs_ext']}"
else:
    CONTIGS_EXT = config["metagenomes"]["contigs_ext"]

# If no contigs given, empty
CONTIGS_NAMES = []

if os.path.isdir(CONTIGS_FOLDER):
    # Get all the files int the contigs folder
    (CONTIGS_NAMES,) = glob_wildcards(
        os.path.join(CONTIGS_FOLDER, "{contigs_files}" + CONTIGS_EXT)
    )

    # Name of the file to check
    CONTIGS_FILES = os.path.join(
        OUTPUT_FOLDER,
        "processed_files",
        "assemblies",
        "merge_contigs.fasta",
    )

    reads = False
else:
    CONTIGS_FILES = os.path.join(
        OUTPUT_FOLDER,
        "processed_files",
        "assemblies",
        "all_contigs.fasta",
    )

    reads=True

# Get fastq folder
FASTQ_FOLDER = config["metagenomes"]["reads_folder"]

# Get fastq separator
FASTQ_SEP = config["metagenomes"]["reads_identifier"]

# Get the samples
(FASTQ_SAMPLE, FASTQ_EXT,) = glob_wildcards(
    os.path.join(
        FASTQ_FOLDER, 
        "{fastq_files}" + FASTQ_SEP + "1" + "{fastq_extension}"
    )
)

# Only get one extension
FASTQ_EXT = FASTQ_EXT[0] if "md5" not in FASTQ_EXT[0] else FASTQ_EXT[0][:-4]

# print(get_final_output(
#             outdir=OUTPUT_FOLDER,
#             reads=reads
#         ))
