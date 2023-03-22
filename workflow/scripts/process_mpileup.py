import sys
from Bio import SeqIO
import pandas as pd

##########################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################

# import pileup output and idxstats output

pileup_dtypes = {
    "name": "string",
    "pos": int,
    "n": "string",
    "number_reads_cov": int,
    "base": "string",
    "quality": "string",
}

pileup = pd.read_table(
    snakemake.input.pileup,
    names=["name", "pos", "n", "number_reads_cov", "base", "quality"],
    dtype=pileup_dtypes,
)

idxstats_dtypes = {
    "name": "string",
    "length": int,
    "no_mapped_reads": int,
    "no_unmapped_reads": int,
}

idxstats = pd.read_table(
    snakemake.input.idxstats,
    names=["name", "length", "no_mapped_reads", "no_unmapped_reads"],
    dtype=idxstats_dtypes,
)

# remove any rows with 0 read coverage from pileup
pileup = pileup[pileup["number_reads_cov"] > 0]

# count number of rows containing a given contigs name in pileup (will be one per base covered) and overwrite pileup df
pileup = pileup.name.value_counts().reset_index(name="covered_bases")
pileup.columns = ["name", "covered_bases"]

# have contingency in the even that the dataframe is empty
if pileup.shape[0] == 0:
    pileup = pd.DataFrame(columns=["name", "covered_bases"])

# drop unnecessary columns from idxstats df
idxstats = idxstats[["name", "length", "no_mapped_reads"]]

# join pileup and idxstats
pileup = pileup.merge(idxstats, on="name").reset_index(drop=True)

# calculate the proportion of bases in the contig with at least 1x coverages
pileup["prop_covered"] = pileup["covered_bases"] / pileup["length"]

# get list of contigs add hit by reads
found = pileup["name"].unique()

# create empty dataframe to house contigs that had no mapped reads
missing = {
    "name": [],
    "covered_bases": [],
    "length": [],
    "no_mapped_reads": [],
    "prop_covered": [],
}

# parse the contig fasta to collect see all contig names and add those missing from the 'found' list to the 'missing' dataframe and set the values of these to 0
parser = SeqIO.parse(snakemake.input.fasta, "fasta")

for record in parser:
    if record.id not in found:
        missing["name"].append(record.id)
        missing["covered_bases"].append(0)
        missing["length"].append(len(record.seq))
        missing["no_mapped_reads"].append(0)
        missing["prop_covered"].append(0.0)
missing = pd.DataFrame.from_dict(missing)

# join the 'missing' dataframe to the one pileup results
pileup = pd.concat([pileup, missing])

# write table to to file
pileup.to_csv(snakemake.output.tsv, sep="\t", index=False, header=True)
