import pandas as pd
import subprocess
import argparse
import os, sys
import subprocess
import pandas as pd

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################


def process_hit_file(blast_path, anicalc_path):

    # import blast hits
    blast_hits = pd.read_csv(
        blast_path,
        sep="\t",
        names=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "qlen",
            "sstart",
            "send",
            "slen",
            "evalue",
            "bitscore",
            "staxids",
            "stitle",
        ],
    )

    # move/drop columns from blast_hits to match required columns for anicalc.py
    blast_hits = blast_hits[
        [
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
            "slen",
        ]
    ]

    # write blast_hits to temp file
    blast_hits.to_csv("temp_blast_hits.tsv", sep="\t", index=False, header=False)

    # run anicalc.py on blast_hits
    subprocess.call(
        [
            "python",
            anicalc_path,
            "-i",
            "temp_blast_hits.tsv",
            "-o",
            "temp_blast_hits.ani",
        ]
    )

    # import blast hits with ANI
    blast_hits = pd.read_csv("temp_blast_hits.ani", sep="\t")

    # filter out hits with ANI < 95% and < 85% coverage of the shortest sequence
    blast_hits = blast_hits[
        (blast_hits["pid"] >= 95)
        & (blast_hits["qcov"] >= 85)
        & (blast_hits["tcov"] >= 85)
    ]

    # delete the 2 temp files
    for f in ["temp_blast_hits.tsv", "temp_blast_hits.ani"]:
        os.remove(f)

    return blast_hits


###########################################################


def process_ictv_hits(ICTV_blast_hits, ICTV_taxonomy):

    # import ICTV taxonomy metadata
    ICTV_taxonomy = pd.read_csv(ICTV_taxonomy, sep="\t")

    # Merge the 2 tables
    ICTV_blast_hits = ICTV_blast_hits.merge(
        ICTV_taxonomy, left_on="tname", right_on="Virus GENBANK accession", how="left"
    )

    # Merge the taxonomy columns into a new column called taxonomy
    ICTV_blast_hits["taxonomy"] = (
        "Viruses;"
        + ICTV_blast_hits["Realm"].astype(str)
        + ";"
        + ICTV_blast_hits["Subrealm"].astype(str)
        + ";"
        + ICTV_blast_hits["Kingdom"].astype(str)
        + ";"
        + ICTV_blast_hits["Subkingdom"].astype(str)
        + ";"
        + ICTV_blast_hits["Phylum"].astype(str)
        + ";"
        + ICTV_blast_hits["Subphylum"].astype(str)
        + ";"
        + ICTV_blast_hits["Class"].astype(str)
        + ";"
        + ICTV_blast_hits["Subclass"].astype(str)
        + ";"
        + ICTV_blast_hits["Order"].astype(str)
        + ";"
        + ICTV_blast_hits["Suborder"].astype(str)
        + ";"
        + ICTV_blast_hits["Family"].astype(str)
        + ";"
        + ICTV_blast_hits["Subfamily"].astype(str)
        + ";"
        + ICTV_blast_hits["Genus"].astype(str)
        + ";"
        + ICTV_blast_hits["Subgenus"].astype(str)
        + ";"
        + ICTV_blast_hits["Species"].astype(str)
    )

    # replace nan with empty string
    ICTV_blast_hits["taxonomy"] = ICTV_blast_hits["taxonomy"].str.replace("nan", "")

    # drop columns from ICTV_blast_hits that arent needed
    ICTV_blast_hits = ICTV_blast_hits[["qname", "taxonomy"]].copy()

    # add column to identify which db the taxon came from
    ICTV_blast_hits["db"] = "ICTV"

    return ICTV_blast_hits


###########################################################


def process_crass_hits(crass_blast_hits):

    # add crass family level taxonomy to hits
    crass_blast_hits[
        "taxonomy"
    ] = "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Crassvirales"

    # drop columns from crass_blast_hits that arent needed
    crass_blast_hits = crass_blast_hits[["qname", "taxonomy"]].copy()

    # add column to identify which db the taxon came from
    crass_blast_hits["db"] = "crass"

    return crass_blast_hits


###########################################################


def process_refseq_hits(refseq_blast_hits, refseq_taxonomy):

    # iterate over each row in refseq_blast_hits and search NCBI for taxonomy using entrez and add the returned taxonomy to the row
    refseq_taxonomy_df = pd.read_table(refseq_taxonomy)

    # Merge the 2 tables
    refseq_blast_hits = refseq_blast_hits.merge(
        refseq_taxonomy_df, left_on="tname", right_on="Representative", how="left"
    )

    refseq_blast_hits = refseq_blast_hits.rename(columns={"Lineage": "taxonomy"})

    # drop columns from refseq_blast_hits that arent needed
    refseq_blast_hits = refseq_blast_hits[["qname", "taxonomy"]].copy()

    # add column to identify which db the taxon came from
    refseq_blast_hits["db"] = "refseq"

    return refseq_blast_hits


###########################################################


def process_imgvr_hits(imgvr_blast_hits, imgvr_taxonomy):

    # fix tname column in imgvr_blast_hits
    imgvr_blast_hits["fixed_tname"] = imgvr_blast_hits["tname"].str.split("|").str[0]

    # import IMG_VR taxonomy metadata
    imgvr_taxonomy = pd.read_csv(
        imgvr_taxonomy,
        sep="\t",
        dtype={"Scaffold_oid": "string", "Host prediction method": "string"},
    )

    # drop columns from imgvr_taxonomy that arent needed
    imgvr_taxonomy = imgvr_taxonomy[["UVIG", "Taxonomic classification"]]

    # Merge the 2 tables
    imgvr_blast_hits = imgvr_blast_hits.merge(
        imgvr_taxonomy, left_on="fixed_tname", right_on="UVIG", how="left"
    )

    # drop columns from imgvr_blast_hits that arent needed
    imgvr_blast_hits = imgvr_blast_hits[["qname", "Taxonomic classification"]]

    # rename columns
    imgvr_blast_hits.columns = ["qname", "taxonomy"]

    # fix taxonomy column in imgvr_blast_hits to be in format "Duplodnaviria;Heunggongvirae;..." rather than "r__Duplodnaviria;k__Heunggongvirae;..."
    imgvr_blast_hits["taxonomy"] = (
        imgvr_blast_hits["taxonomy"]
        .str.replace("r__", "")
        .str.replace("k__", "")
        .str.replace("p__", "")
        .str.replace("c__", "")
        .str.replace("o__", "")
        .str.replace("f__", "")
        .str.replace("g__", "")
        .str.replace("s__", "")
    )

    # append "Viruses;" to the start of each taxonomy
    imgvr_blast_hits["taxonomy"] = "Viruses;" + imgvr_blast_hits["taxonomy"]

    # add column to identify which db the taxon came from
    imgvr_blast_hits["db"] = "IMG_VR"

    return imgvr_blast_hits


###########################################################


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "--ani_script",
        help="The path the the anicalc.py script",
        type=str,
        default=snakemake.params.anicalc,
    )
    parser.add_argument(
        "--outfile",
        help="Path to the output file",
        type=str,
        default=snakemake.output.outfile,
    )
    parser.add_argument(
        "--ictv_blastout",
        help="The blast outfmt6 file against ICTV",
        type=str,
        default=snakemake.input.ictv_blasts,
    )
    parser.add_argument(
        "--refseq_blastout",
        help="The blast outfmt6 file against RefSeq Viral",
        type=str,
        default=snakemake.input.refseq_blasts,
    )
    parser.add_argument(
        "--crass_blastout",
        help="The blast outfmt6 file against Crassphage database",
        type=str,
        default=snakemake.input.crass_blasts,
    )
    parser.add_argument(
        "--img_vr_blastout",
        help="The blast outfmt6 file against IMG_VR",
        type=str,
        default=snakemake.input.img_blasts,
    )
    parser.add_argument(
        "--ictv_db",
        help="Path to the ICTV metadata tabel with taxonomy",
        type=str,
        default=snakemake.params.ictv_db,
    )
    parser.add_argument(
        "--refseq_db",
        help="Path to the RefSeq Viral metadata tabel with taxonomy",
        type=str,
        default=snakemake.params.refseq_db,
    )
    parser.add_argument(
        "--imgvr_db",
        help="Path to the IMG_VR metadata tabel with taxonomy",
        type=str,
        default=snakemake.params.img_db,
    )

    # Parse arguments
    args = parser.parse_args()

    return args


###########################################################

args = parse_args()

# get species level hit ICTV contigs
ICTV_blast_hits = process_hit_file(args.ictv_blastout, args.ani_script)

# get taxa for the good hit ICTV contigs
ICTV_blast_hits = process_ictv_hits(ICTV_blast_hits, args.ictv_db)

# get species level hit refseq contigs
refseq_blast_hits = process_hit_file(args.refseq_blastout, args.ani_script)

# get taxa for the good hit refseq contigs
refseq_blast_hits = process_refseq_hits(refseq_blast_hits, args.refseq_db)

# get species level hit crass contigs
crass_blast_hits = process_hit_file(args.crass_blastout, args.ani_script)

# get taxa for the good hit crass contigs
crass_blast_hits = process_crass_hits(crass_blast_hits)

# get species level hit img_vr contigs
imgvr_blast_hits = process_hit_file(args.img_vr_blastout, args.ani_script)

# get taxa for the good hit img_vr contigs
imgvr_blast_hits = process_imgvr_hits(imgvr_blast_hits, args.imgvr_db)

# combine all the hits into one table
all_hits = pd.concat(
    [ICTV_blast_hits, refseq_blast_hits, crass_blast_hits, imgvr_blast_hits]
)

# identify duplicate qnames and store them in a list
duplicates = all_hits[all_hits.duplicated(["qname"], keep=False)]["qname"].unique().tolist()


# cycle through the list of duplicates and collect hits to dbs in the following order of preference ictv, refseq, img_vr, crass

best_taxa = []

for i in duplicates:
    if len(all_hits[(all_hits["qname"] == i) & (all_hits["db"] == "ICTV")]) > 0:
        best_taxa.append(
            all_hits[(all_hits["qname"] == i) & (all_hits["db"] == "ICTV")]
        )
    elif len(all_hits[(all_hits["qname"] == i) & (all_hits["db"] == "refseq")]) > 0:
        best_taxa.append(
            all_hits[(all_hits["qname"] == i) & (all_hits["db"] == "refseq")]
        )
    elif len(all_hits[(all_hits["qname"] == i) & (all_hits["db"] == "IMG_VR")]) > 0:
        best_taxa.append(
            all_hits[(all_hits["qname"] == i) & (all_hits["db"] == "IMG_VR")]
        )
    elif len(all_hits[(all_hits["qname"] == i) & (all_hits["db"] == "crass")]) > 0:
        best_taxa.append(
            all_hits[(all_hits["qname"] == i) & (all_hits["db"] == "crass")]
        )

best_taxa = pd.concat(best_taxa).reset_index(drop=True)

# keep only the first of perfect duplicates
best_taxa = best_taxa.drop_duplicates(["qname", "db"], keep="first").reset_index(drop=True)

# flag if there are still any remaining duplicates
if len(best_taxa[best_taxa.duplicated(["qname"], keep=False)]["qname"].unique().tolist()) > 0:
    print(
        "There are still duplicates in the best_taxa table fix script as likely conflicting taxonomy was found for the same contig"
    )
    exit(1)

# concat best_taxa with the unique values in the all_hits table
all_hits = all_hits[~all_hits["qname"].isin(duplicates)].reset_index(drop=True)
all_hits = pd.concat([all_hits, best_taxa], ignore_index=True)

# write all_hits to file
all_hits.to_csv(args.outfile, sep="\t", index=False)
