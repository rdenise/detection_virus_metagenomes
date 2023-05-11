import argparse
import pandas as pd
import sys

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

def parse_arguments():
    # create the argument parser
    parser = argparse.ArgumentParser(
        description="Classify bacphlip entries as virulent or temperate based on a cutoff value."
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        help="the cutoff value for classification",
        default=snakemake.params.cutoff,
    )
    parser.add_argument(
        "--input_file",
        help="the path to the input file",
        default=snakemake.input.bacphlip,
    )
    parser.add_argument(
        "--output_file",
        help="the path to the output file",
        default=snakemake.output.modify,
    )

    # parse the arguments
    args = parser.parse_args()

    return args

###########################################################
###########################################################


def main(args):

    # import bacphlip dataframe
    df = pd.read_csv(args.input_file, sep="\t")

    # rename columns
    df.columns = ["contig", "virulent", "temperate"]

    # create a result column that just says which ones are virulent vs temperate based on the argv[1] cutoff
    df["result"] = df.apply(
        lambda x: "virulent"
        if x["virulent"] >= args.cutoff and x["temperate"] < args.cutoff
        else (
            "temperate"
            if x["temperate"] >= args.cutoff and x["virulent"] < args.cutoff
            else "undetermined"
        ),
        axis=1,
    )

    # write new table to file
    df.to_csv(args.output_file, sep="\t", index=False)

###########################################################
###########################################################


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
