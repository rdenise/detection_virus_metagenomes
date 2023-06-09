from Bio import SeqIO
import argparse
import os
import subprocess
import io, sys
# import multiprocessing as mp
import shutil
# import pandas as pd
from pathlib import Path
import time
from concurrent.futures import ProcessPoolExecutor

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################


def create_folder(mypath):

    """
    Created the folder that I need to store my result if it doesn't exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return


###########################################################


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "--in_fasta_file", help="input_fasta", type=str, default=snakemake.input.contig
    )
    parser.add_argument(
        "--outfile", help="output_folder", type=str, default=snakemake.output.blast_out
    )
    parser.add_argument(
        "--tmp_dir", help="tmp folder", type=str, default=snakemake.params.tmp_output
    )
    parser.add_argument(
        "--evalue", help="e_value", type=str, default=snakemake.params.evalue
    )
    parser.add_argument(
        "--outfmt", help="output blast format", type=str, default=snakemake.params.outfmt
    )
    parser.add_argument(
        "--max_length", help="max length of the sequence", type=str, default=snakemake.params.max_len
    )    
    parser.add_argument(
        "--database",
        help="The name of the database or empty if remote blast",
        type=str,
        default=snakemake.params.database,
    )
    parser.add_argument(
        "--contigs_per_file",
        help="The number of contigs in each slice",
        type=int,
        default=1050,
    )
    parser.add_argument("--job_number", type=int, default=snakemake.threads)
    parser.add_argument(
        "--options_blast", help="other blast option", type=str, default=snakemake.params.options_blast
    )
    # Parse arguments
    args = parser.parse_args()

    return args


###########################################################


def main(args):

    # split the input file
    record_iter = open(args.in_fasta_file)
    files_to_run = []

    for i, batch in enumerate(batch_iterator(record_iter, args.contigs_per_file, args.max_length)):
        try:
            group_name = f"group_{i+1}"
            dir_name = f"{args.tmp_dir}/{group_name}"
            filename = f"{dir_name}/{group_name}.fasta"

            create_folder(f"{args.tmp_dir}/{group_name}")

            with open(filename, "w") as handle:
                # count = SeqIO.write(batch, handle, "fasta")
                handle.write(batch[0])
                print(f"Wrote {batch[1]} records to {filename}")
                files_to_run.append((group_name, dir_name, filename))
        except IOError as ioerror:
            print(ioerror)

    if args.database:
        if not (
            Path(args.database + ".00.nsq").exists()
            or Path(args.database + ".nsq").exists()
        ):
            cmd_str = f"makeblastdb -dbtype nucl -in '{args.database}'"
            stdout, stderr = execute(cmd_str)
            print(
                f"----makeblastdb - stdout----\n{stdout}\n----makeblastdb - stderr----\n{stderr}\n"
            )

    record_iter.close()

    # pool = mp.Pool(args.job_number)
    # results = pool.map(run_job, files_to_run)
    # pool.close()

    with ProcessPoolExecutor(max_workers=args.job_number) as executor:
        results = list(executor.map(run_job, files_to_run))

    # results = Parallel(n_jobs=args.job_number)(delayed(run_job)(file) for file in files_to_run)

    # df = pd.concat(results)
    # df.to_csv(args.outfile, sep="\t", index=False, header=False)

    with open(args.outfile, "wb") as outfile:
        for result in results:
            with open(os.path.join(result), 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

    print(f"{bcolors.OKBLUE} ------ DONE! ----------- {bcolors.ENDC}")
    shutil.rmtree(args.tmp_dir)


###########################################################


def run_job(group_tuple):
    if args.database:
        blast_remote = ""
        blast_database = args.database
    else:
        blast_remote = "-task blastn -remote"
        blast_database = "nt"

    # -out {group_tuple[1]}/blast-output.txt 

    job_str = (
        f"blastn -query {group_tuple[2]} -out {group_tuple[1]}/blast-output.txt "
        f"-db '{blast_database}' -evalue {args.evalue} {blast_remote} "
        f"-outfmt '{args.outfmt}' {args.options_blast}"
    )

    file_name = f"{group_tuple[1]}/blast-output.txt"

    if os.path.isfile(file_name):
        stdout = f"File {file_name} already exists"
        stderr = ""
    else:
        stdout, stderr = execute(job_str)
    # stdout, stderr = execute(job_str)

    print(f"----BLASTn - stderr----\n{stderr}\n")
    # print(f"----BLASTn - stdout----\n{stdout}\n----BLASTn - stderr----\n{stderr}\n")

    # if os.path.isfile(file_name) and os.path.getsize(
    #     file_name
    # ):
    #     df = pd.read_csv(file_name, sep="\t", header=None)
    # else:
    #     df = pd.DataFrame()

    return file_name


###########################################################


def execute(command):
    print(f"Executing {command}")
    process = subprocess.Popen(
        command,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
    )
    stdout, stderr = process.communicate()
    return stdout, stderr


###########################################################


def batch_iterator(iterator, batch_size, max_len):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        index = 0
        while index < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            # if max_len and len(entry.seq) <= max_len:
            elif entry.startswith(">"):
                index += 1
            
            batch.append(entry)
        if batch:
            yield "".join(batch[:-1]), index
            batch = [entry]


###########################################################


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


###########################################################

if __name__ == "__main__":
    args = parse_args()
    main(args)
