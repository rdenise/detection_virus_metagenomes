import sys
from pathlib import Path
import pandas as pd

##########################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

##########################################################################

# import total number of reads in the samples
total_reads_dtypes={
    "sample":"string", 
    "total_reads":int,
}

total_reads_count_df = pd.read_table(snakemake.input.total_reads, dtype=total_reads_dtypes)
total_reads_count_dict = total_reads_count_df.set_index("sample").read_no.to_dict()

# import all mpileup files and use to create new pandas dataframe
mpileup_df = []

for mpileup_cov in snakemake.input.mpileup_cov:
    sub_df = pd.read_table(mpileup_cov)
    sub_df["sample"] = Path(mpileup_cov).stem.split("_mpileup")[0]
    mpileup_df.append(sub_df)

mpileup_df = pd.concat(mpileup_df, ignore_index=True)

# process mpileup df to set coverages to 0 for any contig that less than 75% of its length has been covered by at least 1x coverage
mpileup_df["width_adjusted_counts"] = mpileup_df["no_mapped_reads"]
mpileup_df.loc[mpileup_df.prop_covered < 0.75, "width_adjusted_counts"] = 0

# normalise read counts by contig length (i.e. counts/base)
mpileup_df["length_normalised_counts"] = mpileup_df["width_adjusted_counts"] / mpileup_df["length"]

# total-sum normalise each read set by dividing the normalised reads by the total number of quality filtered reads in the sample (could be multiprocessed)
# chris = total_reads
mpileup_df["total_reads"] = mpileup_df["sample"].map(total_reads_count_dict)
mpileup_df["total_sum_normalised"] = mpileup_df.length_normalised_counts / mpileup_df.total_reads

# replace any nan values in the total-sum normalised column with 0s
mpileup_df["total_sum_normalised"] = mpileup_df["total_sum_normalised"].fillna(0)

# output the mpileup_df to a file (chris = tss_counts)
mpileup_df.to_csv(snakemake.output.long_format, sep = "\t", index = False)

# output total-sum_normalised otu table
mpileup_df.pivot(index="name", columns="sample", values= "total_sum_normalised").to_csv(snakemake.output.otu, sep = "\t", index = False)



