from Bio import SeqIO
import pandas as pd
import sys
import os

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

# Dataframe that contains all the informations about
genomad_contigs = pd.read_table(snakemake.input.virus_summary).seq_name.unique().tolist()

all_contigs = genomad_contigs

with open(snakemake.output.fasta, "wt") as w_file:
    parser = SeqIO.parse(snakemake.input.fasta, "fasta")

    for contig in parser:
        if contig.id not in all_contigs:
            SeqIO.write(contig, w_file, "fasta")

###########################################################
###########################################################
