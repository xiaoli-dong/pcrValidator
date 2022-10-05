#!/usr/bin/env python3

# pcrValidator- A tool for assessing PCR assay through the evaluation of mismatches with publicly available sequences

# Xiaoli Dong
# Public Health Laboratory (ProvLab) - Alberta Health Services
# Xiaoli.Dong@albertaprecisionlabs.ca

import subprocess
import os
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import AlignIO

import distance

__version__ = "0.1.0"
__author__ = "Xiaoli Dong"


def get_parser():

    # Disable default help
    parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="For example:\n python pcrValidator.py -a assay.csv -g mpx.fasta -o results -p mpx -c 0",
    )
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    # Add back help
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit",
    )

    required.add_argument(
        "-a", type=str, required=True, help="path to assay details in CSV file"
    )
    required.add_argument(
        "-g", type=str, required=True, help="path to target genomes in FASTA file"
    )

    # prog give the program basename without path but __file__ gives the full path
    version = parser.prog + " " + __version__
    optional.add_argument("--version", "-v", action="version", version=version)

    optional.add_argument(
        "--outdir",
        "-o",
        type=str,
        default=".",
        help="path to output directory",
    )
    optional.add_argument(
        "--prefix",
        "-p",
        type=str,
        default="output",
        help="output file name prefix",
    )
    optional.add_argument(
        "--cutoff",
        "-c",
        type=float,
        default=0,
        help="minimum prevalence (%%) of primer site variants reported in reports, min > 0, max < 100",
    )
    optional.add_argument(
        "-e",
        type=int,
        default=45,
        help="minimum primer Tm (degrees C)",
    )
    optional.add_argument(
        "-E",
        type=int,
        default=45,
        help="minimum probe Tm (degrees C)",
    )
    optional.add_argument(
        "-t",
        type=float,
        default=0.9,
        help="primer strand concentration (in uM)",
    )
    optional.add_argument(
        "-T",
        type=float,
        default=0.25,
        help="Probe strand concentration (in uM)",
    )

    return parser


def check_genomes_file(path_to_file):
    """Checks that file exists, is not empty, and headers are unique and hashable."""
    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f"\nERROR: Genomes FASTA file {path_to_file} does not exist.\n")
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f"\nERROR: Genomes FASTA file {path_to_file} is empty.\n")
        exit(1)
    with open(path_to_file, "r") as input_file:
        headers = []
        for line in input_file:
            if line[0] == ">":
                headers.append(line.strip().lstrip(">"))
        if len(headers) != len(set(headers)):
            print("\nERROR: FASTA headers in input genomes file must be unique.\n")
            exit(1)


def read_assay_file(path_to_file):
    """Read input assay file and parse each list into a tuple:
    (assay_name,fwd_oligo_name,fwd_oligo_seq,rev_oligo_name,rev_oligo_seq,probe_oligo_name,probe_oligo_seq)
    Return a tuple of these tuples."""
    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(f"\nERROR: Assay CSV file {path_to_file} does not exist.\n")
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(f"\nERROR: Assay CSV file {path_to_file} is empty.\n")
        exit(1)
    assays = []
    with open(path_to_file, "r") as input_file:
        for line in input_file:
            line = line.rstrip().split(",")
            # Check if line has the accepted number of comma-separated values
            if len(line) not in [5, 7]:
                print("\nERROR: Line in assay file not properly formatted:")
                print(",".join(line))
                print("\nAssay file lines must be a comma-separated list:")
                print(
                    "assay_name,fwd_primer_name,fwd_primer_seq,rev_primer_name,rev_primer_seq,probe_name,probe_seq"
                )
                exit(1)
            else:
                assay_name = line[0]
                fwd_primer_name = line[1]
                fwd_primer_seq = line[2]
                rev_primer_name = line[3]
                rev_primer_seq = line[4]
                if len(line) == 7:
                    probe_name = line[5]
                    probe_seq = line[6]
                else:
                    probe_name = ""
                    probe_seq = ""
                # Check if any assay or primer names are empty
                names = [assay_name, fwd_primer_name, rev_primer_name]
                names += [probe_name] if probe_seq != "" else []
                if any([name == "" for name in names]):
                    print("\nERROR: Assay and oligo names cannot be empty:")
                    print(",".join(line))
                    print()
                    exit(1)
                # Check if oligo seqs contain invalid bases
                bases = "ATGCWSMKRYBVDHN"
                for oligo, seq in zip(
                    ["fwd primer", "rev primer", "probe"],
                    [fwd_primer_seq, rev_primer_seq, probe_seq],
                ):
                    if len(set(seq.upper()) - set(bases)) > 0:
                        print(f"\nERROR: Oligo seq for {oligo} contains invalid bases:")
                        print(seq)
                        print()
                        exit(1)
                # Create assay info tuple and append to assays list
                assay = (
                    assay_name,
                    fwd_primer_name,
                    fwd_primer_seq,
                    rev_primer_name,
                    rev_primer_seq,
                    probe_name,
                    probe_seq,
                )
                assays.append(assay)
    # Check that names used for assays and oligos are unique
    for name_index in [0, 1, 3, 5]:
        names = [assay[name_index] for assay in assays if assay[name_index] != ""]
        if len(names) != len(set(names)):
            print("\nERROR: Assay and oligo names must be unique!\n")
            exit()
    return tuple(assays)


def run_TNTBLAST(
    assay_details,
    path_to_genomes,
    outdir,
    prefix,
    melting_temp_primer,
    melting_temp_probe,
    primer_molarity,
    probe_molarity,
):
    """Run TNTBLAST with provided assay details."""
    # Create TNTBLAST assay file from assay details (TSV file of assay name and oligo seqs)
    assay_name = assay_details[0]
    fwd_primer_seq = assay_details[2]
    rev_primer_seq = assay_details[4]
    probe_seq = assay_details[6]
    print(f"Running TNTBLAST for {assay_name} against {path_to_genomes}...")
    print(f"Fwd primer seq: {fwd_primer_seq}\nRev primer seq: {rev_primer_seq}")
    if probe_seq != "":
        print(f"Probe seq: {probe_seq}")
    assay = [assay_name, fwd_primer_seq, rev_primer_seq]
    assay += [probe_seq] if probe_seq != "" else []
    assay = "\t".join(assay)
    path_to_tntblast_assay = os.path.join(outdir, assay_name + "_tntblast_assay.tsv")
    with open(path_to_tntblast_assay, "w") as output_file:
        output_file.write(assay + "\n")
    # Create terminal command for TNTBLAST and run
    path_to_tntblast_txt_output = os.path.join(
        outdir, prefix + "_" + assay_name + "_amplicon.fasta"
    )
    terminal_command = (
        f"tntblast -i {path_to_tntblast_assay} -d {path_to_genomes} -o {path_to_tntblast_txt_output}"
        f" -e {melting_temp_primer} -E {melting_temp_probe} -t {primer_molarity / 1000000} -T {probe_molarity / 1000000}"
        f" --best-match -m 1 -v F --single-primer-pcr F"
    )
    print(terminal_command)
    completed_process = subprocess.run(
        terminal_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True,
    )
    if completed_process.returncode != 0:
        print(
            f"\nERROR: TNTBLAST terminated with errors.\nTNTBLAST error code: {completed_process.returncodes}\n"
        )
        exit(1)
    # GARBAGE COLLECTION (remove TNTBLAST assay TSV file)
    os.remove(path_to_tntblast_assay)
    print()
    return path_to_tntblast_txt_output


def count_targets(path_to_genomes):
    """Counts the number of target sequences in the provided genomes FASTA file.
    Returns an int of the count."""
    with open(path_to_genomes, "r") as input_file:
        seq_counter = 0
        for line in input_file:
            if line[0] == ">":
                seq_counter += 1
    return int(seq_counter)


def get_targets(path_to_genomes):
    """Reads all sequence in the provided genomes FASTA file and returns them in
    a dict where keys are FASTA headers and values are the nucleotide sequences."""
    with open(path_to_genomes, "r") as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == ">":
                header = line.strip().lstrip(">")
                seqs[header] = ""
            else:
                seqs[header] += line.strip()
    return seqs


def run_msa(assay_details, prefix, outdir):

    # Write primer probe to a file and convert reverse primer to be on the plus strand

    assay_name = assay_details[0]
    fwd_primer_name = assay_details[1]
    fwd_primer_seq = assay_details[2]
    rev_primer_name = assay_details[3]
    rev_primer_seq = Seq(assay_details[4]).reverse_complement()
    probe_name = assay_details[5]
    probe_seq = assay_details[6]

    if probe_seq != "":
        print(f"Probe seq: {probe_seq}")

    assay = [[fwd_primer_name, fwd_primer_seq], [rev_primer_name, rev_primer_seq]]
    # assay += [probe_name, probe_seq] if probe_seq != '' else []
    if probe_seq != "":
        assay.append([probe_name, probe_seq])

    df = pd.DataFrame(assay, columns=["Name", "seq"])
    df["Name"] = ">" + df["Name"]

    # the orientiation of the sequences is on the plus strand
    path_to_assay_seq = os.path.join(
        outdir, prefix + "_" + assay_name + "_msa_assay.fasta"
    )
    df.to_csv(path_to_assay_seq, sep="\n", index=False, header=False)

    """Run mafft with provided assay details."""
    path_to_amplicon = os.path.join(
        outdir, prefix + "_" + assay_name + "_amplicon.fasta"
    )
    print(
        f"Running mafft to produce amplicon multiple sequence alignment for {path_to_amplicon}..."
    )
    # Create terminal command for TNTBLAST and run
    path_to_mafft_output = os.path.join(
        outdir, prefix + "_" + assay_name + "_amplicon_mafft.fasta"
    )
    terminal_command = f"cat {path_to_assay_seq} {path_to_amplicon} | mafft --auto --thread -1 - > {path_to_mafft_output}"
    print(terminal_command)
    completed_process = subprocess.run(
        terminal_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True,
    )
    if completed_process.returncode != 0:
        print(
            f"\nERROR: mafft terminated with errors.\nmafft error code: {completed_process.returncodes}\n"
        )
        exit(1)
    # GARBAGE COLLECTION (remove TNTBLAST assay TSV file)
    # os.remove(path_to_tntblast_assay)
    print()


def parse_msa_output(assay_details, prefix, outdir, path_to_genomes):

    assay_name = assay_details[0]
    fwd_primer_name = assay_details[1]
    rev_primer_name = assay_details[3]
    probe_name = assay_details[5]

    ################################
    # path_to_tntblast_fasta_output = os.path.join(
    #     outdir, assay_name + "_amplicon.fasta"
    # )
    path_to_tntblast_fasta_output = os.path.join(
        outdir, prefix + "_" + assay_name + "_amplicon.fasta"
    )
    # Count targets detected by each assay

    # targets_detected = tntblast_results["assay_name"].value_counts()
    total_input_seqs = count_targets(path_to_genomes)
    amplicon_count = count_targets(path_to_tntblast_fasta_output)
    #################################

    path_to_mafft_output = os.path.join(
        # outdir, prefix + "_" + assay_name + "_mafft_output.fasta"
        outdir,
        prefix + "_" + assay_name + "_amplicon_mafft.fasta",
    )
    align = AlignIO.read(path_to_mafft_output, "fasta")

    # print(type(align))
    # print(align)
    fwd_coord_found = False
    rev_coord_found = False
    probe_coord_found = False

    fwd_primer_aligned_str = ""
    rev_primer_aligned_str = ""
    probe_aligned_str = ""

    msa_results = pd.DataFrame()
    fwd_df = pd.DataFrame()
    rev_df = pd.DataFrame()
    probe_df = pd.DataFrame()

    for record in align:

        if record.id == fwd_primer_name:
            fwd_primer_alignment = str(record.seq)
            left_stripped = fwd_primer_alignment.lstrip("-")
            right_stripped = left_stripped.rstrip("-")
            fwd_start = len(fwd_primer_alignment) - len(left_stripped)
            fwd_end = len(right_stripped) + fwd_start
            fwd_primer_aligned_str = right_stripped
            fwd_coord_found = True
            fwd_primer_site_list = []
            target_list = []

            for record in align[:, fwd_start:fwd_end]:
                # print(str(fwd_record.seq))
                if "n" in str(record.seq):
                    continue
                elif (
                    (record.id == fwd_primer_name)
                    or (record.id == rev_primer_name)
                    or (record.id == probe_name)
                ):
                    continue
                fwd_primer_site_list.append(str(record.seq))
                target_list.append(record.id)

            fwd_df = pd.DataFrame(
                list(zip(target_list, fwd_primer_site_list)),
                columns=["target", "fwd_site"],
            )

            fwd_df["fwd_errors"] = fwd_df.apply(
                lambda x: distance.hamming(x["fwd_site"], fwd_primer_aligned_str),
                axis=1,
            )

        elif record.id == rev_primer_name:
            rev_primer_alignment = str(record.seq)
            left_stripped = rev_primer_alignment.lstrip("-")
            right_stripped = left_stripped.rstrip("-")
            rev_start = len(rev_primer_alignment) - len(left_stripped)
            rev_end = len(right_stripped) + rev_start
            rev_primer_aligned_str = right_stripped
            rev_coord_found = True
            rev_primer_site_list = []
            target_list = []
            for record in align[:, rev_start:rev_end]:
                # print(str(fwd_record.seq))
                if "n" in str(record.seq):
                    continue
                elif (
                    (record.id == fwd_primer_name)
                    or (record.id == rev_primer_name)
                    or (record.id == probe_name)
                ):
                    continue
                rev_primer_site_list.append(str(record.seq))
                target_list.append(record.id)

            rev_df = pd.DataFrame(
                list(zip(target_list, rev_primer_site_list)),
                columns=["target", "rev_site"],
            )

            rev_df["rev_errors"] = rev_df.apply(
                lambda x: distance.hamming(x["rev_site"], rev_primer_aligned_str),
                axis=1,
            )

            # print(rev_df)

        elif probe_name != "" and record.id == probe_name:
            probe_alignment = str(record.seq)
            left_stripped = probe_alignment.lstrip("-")
            right_stripped = left_stripped.rstrip("-")
            probe_start = len(probe_alignment) - len(left_stripped)
            probe_end = len(right_stripped) + probe_start
            probe_aligned_str = right_stripped
            probe_coord_found = True
            probe_site_list = []
            target_list = []
            for record in align[:, probe_start:probe_end]:
                # print(str(fwd_record.seq))
                if "n" in str(record.seq):
                    continue
                elif (
                    (record.id == fwd_primer_name)
                    or (record.id == rev_primer_name)
                    or (record.id == probe_name)
                ):
                    continue
                probe_site_list.append(str(record.seq))
                target_list.append(record.id)

            probe_df = pd.DataFrame(
                list(zip(target_list, probe_site_list)),
                columns=["target", "probe_site"],
            )

            probe_df["probe_errors"] = probe_df.apply(
                lambda x: distance.hamming(x["probe_site"], probe_aligned_str),
                axis=1,
            )

            # print(probe_df)

        ################## producing report , the assay has no probe

        if assay_details[5] == "" and fwd_coord_found and rev_coord_found:
            # msa_results = pd.merge(fwd_df_mask, rev_df_mask, on="target")
            msa_results = pd.merge(fwd_df, rev_df, on="target")
            msa_results["total_errors"] = (
                msa_results["fwd_errors"] + msa_results["rev_errors"]
            )
            msa_results = msa_results[
                (msa_results.target != fwd_primer_name)
                & (msa_results.target != rev_primer_name)
            ]

            msa_results.to_csv(
                os.path.join(
                    outdir,
                    prefix + "_" + assay_name + "_pcr.tsv",
                ),
                sep="\t",
                index=False,
                # index_label="variant_name",
            )

            # print(msa_results)
            break
        elif (
            assay_details[5] != ""
            and fwd_coord_found
            and rev_coord_found
            and probe_coord_found
        ):
            msa_results = pd.merge(fwd_df, probe_df, on="target")
            msa_results = pd.merge(msa_results, rev_df, on="target")

            msa_results["total_errors"] = (
                msa_results["fwd_errors"]
                + msa_results["rev_errors"]
                + msa_results["probe_errors"]
            )

            msa_results.to_csv(
                os.path.join(
                    outdir,
                    prefix + "_" + assay_name + "_pcr.tsv",
                ),
                sep="\t",
                index=False,
                # index_label="variant_name",
            )

            # print(msa_results)
            break

    msa_results["assay_name"] = assay_name
    msa_results["total_input_seqs"] = total_input_seqs
    msa_results["amplicon_count"] = amplicon_count

    return (msa_results,)


def write_assay_report(assay_details, prefix, outdir, msa_results, threshold):

    assay_name = assay_details[0]

    # Count targets detected by each assay
    # targets_detected = msa_results["assay_name"].value_counts()
    # total_input_seqs = count_targets(path_to_genomes)
    # amplicon_count = targets_detected[0]
    total_input_seqs = msa_results["total_input_seqs"].iloc[0]
    amplicon_count = msa_results["amplicon_count"].iloc[0]
    # taqman assay with probe
    if assay_details[5] != "":
        all_cols = [
            "fwd_site",
            "fwd_errors",
            "probe_site",
            "probe_errors",
            "rev_site",
            "rev_errors",
            "total_errors",
        ]
        all_summary = (
            msa_results.groupby(all_cols).size().reset_index(name="all_site_count")
        )

        all_summary["total_input_seq"] = total_input_seqs
        all_summary["total_amplicon_count"] = amplicon_count
        all_summary["amplicon_in_input_pct"] = round(
            amplicon_count * 100 / total_input_seqs, 2
        )
        all_summary["all_site_in_amplicon_pct"] = round(
            all_summary["all_site_count"] * 100 / amplicon_count, 2
        )

        all_summary.astype(
            {
                "total_errors": "int",
                "all_site_count": "int",
                "total_amplicon_count": "int",
                "total_input_seq": "int",
            }
        )
        all_summary = all_summary[
            all_summary["all_site_in_amplicon_pct"] >= float(threshold)
        ]
        all_summary.sort_values(
            ["all_site_count"], axis=0, ascending=[False], inplace=True
        )
        all_summary.to_csv(
            os.path.join(
                outdir,
                prefix + "_" + assay_name + "_assay_report.tsv",
            ),
            sep="\t",
            index=False,
            # index_label="variant_id",
        )
    else:
        # print(msa_results)
        all_cols = ["fwd_site", "fwd_errors", "rev_site", "rev_errors", "total_errors"]
        all_summary = (
            msa_results.groupby(all_cols).size().reset_index(name="all_site_count")
        )

        all_summary["total_input_seq"] = total_input_seqs
        all_summary["total_amplicon_count"] = amplicon_count
        all_summary["amplicon_in_input_pct"] = round(
            amplicon_count * 100 / total_input_seqs, 2
        )
        all_summary["all_site_in_amplicon_pct"] = round(
            all_summary["all_site_count"] * 100 / amplicon_count, 2
        )

        all_summary.astype(
            {
                "total_errors": "int",
                "all_site_count": "int",
                "total_amplicon_count": "int",
                "total_input_seq": "int",
            }
        )
        all_summary = all_summary[
            all_summary["all_site_in_amplicon_pct"] >= float(threshold)
        ]
        all_summary.sort_values(
            ["all_site_count"], axis=0, ascending=[False], inplace=True
        )
        all_summary.to_csv(
            os.path.join(
                outdir,
                prefix + "_" + assay_name + "_assay_report.tsv",
            ),
            sep="\t",
            index=False,
            # index_label="variant_id",
        )


def write_variants_report(assay_details, prefix, outdir, msa_results, threshold):

    assay_name = assay_details[0]

    # Count targets detected by each assay
    targets_detected = msa_results["assay_name"].value_counts()
    amplicon_count = targets_detected[0]
    total_input_seqs = msa_results["total_input_seqs"].iloc[0]
    amplicon_count = msa_results["amplicon_count"].iloc[0]
    assay_name = msa_results["assay_name"].iloc[0]

    fwd_summary = (
        msa_results.groupby(["fwd_site", "fwd_errors"])
        .size()
        .reset_index(name="fwd_site_count")
    )

    # fwd_summary["total_input_seq"] = total_input_seqs
    fwd_summary["total_amplicon_count"] = amplicon_count
    # fwd_summary["amplicon_in_input_pct"] = round(
    #     targets_detected[0] * 100 / total_input_seqs, 2
    # )
    fwd_summary["fwd_in_amplicon_pct"] = round(
        fwd_summary["fwd_site_count"] * 100 / amplicon_count, 2
    )

    fwd_summary.astype(
        {
            "fwd_errors": "int",
            "fwd_site_count": "int",
            "total_amplicon_count": "int",
            # "total_input_seq": "int",
        }
    )
    fwd_summary = fwd_summary[fwd_summary["fwd_in_amplicon_pct"] >= float(threshold)]
    # sort data frame
    fwd_summary.sort_values(["fwd_site_count"], axis=0, ascending=[False], inplace=True)
    fwd_summary.to_csv(
        os.path.join(
            outdir,
            prefix + "_" + assay_name + "_fwd_variants.tsv",
        ),
        sep="\t",
        index=False,
        # index_label="variant_id",
    )
    # print(final_fwd_summary)

    rev_summary = (
        msa_results.groupby(["rev_site", "rev_errors"])
        .size()
        .reset_index(name="rev_site_count")
    )

    # rev_summary["total_input_seq"] = total_input_seqs
    rev_summary["total_amplicon_count"] = amplicon_count
    # rev_summary["amplicon_in_input_pct"] = round(
    #     targets_detected[0] * 100 / total_input_seqs, 2
    # )
    rev_summary["rev_in_amplicon_pct"] = round(
        rev_summary["rev_site_count"] * 100 / amplicon_count, 2
    )

    rev_summary.astype(
        {
            "rev_errors": "int",
            "rev_site_count": "int",
            "total_amplicon_count": "int",
            # "total_input_seq": "int",
        }
    )
    rev_summary = rev_summary[rev_summary["rev_in_amplicon_pct"] >= float(threshold)]
    rev_summary.sort_values(["rev_site_count"], axis=0, ascending=[False], inplace=True)
    rev_summary.to_csv(
        os.path.join(
            outdir,
            prefix + "_" + assay_name + "_rev_variants.tsv",
        ),
        sep="\t",
        index=False,
        # index_label="variant_id",
    )
    # print(final_rev_summary)

    if assay_details[5] != "":
        probe_summary = (
            msa_results.groupby(["probe_site", "probe_errors"])
            .size()
            .reset_index(name="probe_site_count")
        )

        # probe_summary["total_input_seq"] = total_input_seqs
        probe_summary["total_amplicon_count"] = amplicon_count
        # probe_summary["amplicon_in_input_pct"] = round(
        #     targets_detected[0] * 100 / total_input_seqs, 2
        # )
        probe_summary["probe_in_amplicon_pct"] = round(
            probe_summary["probe_site_count"] * 100 / amplicon_count, 2
        )
        probe_summary.astype(
            {
                "probe_errors": "int",
                "probe_site_count": "int",
                "total_amplicon_count": "int",
                # "total_input_seq": "int",
            }
        )
        probe_summary = probe_summary[
            probe_summary["probe_in_amplicon_pct"] >= float(threshold)
        ]
        probe_summary.sort_values(
            ["probe_site_count"], axis=0, ascending=[False], inplace=True
        )
        probe_summary.to_csv(
            os.path.join(
                outdir,
                prefix + "_" + assay_name + "_probe_variants.tsv",
            ),
            sep="\t",
            index=False,
            # index_label="variant_id",
        )
        # print(probe_summary)


def main():

    # Parse command line arguments
    # __file__ has full path associate with it,
    parser = get_parser()
    args = parser.parse_args()
    print(f"\n {parser.prog} v{__version__}\n")

    # check output dir exists otherwise create
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    if args.cutoff < 0 or args.cutoff > 100:
        print(
            f"\nERROR: minimum prevalence (%) of primer site variants should in the range [0,100]\n"
        )
        exit(1)
    if args.t < 0:
        print(f"\nmolar concentration of primer oligos (MOL) must be > 0\n")
        exit(1)
    if args.T < 0:
        print(f"\nmolar concentration of probe oligos (MOL) must be > 0\n")
        exit(1)

    # Check input genomes file
    check_genomes_file(args.g)

    # Parse assay file to get oligo names and seqs
    assays = read_assay_file(args.a)

    # Create empty dataframe for tabulated tntblast results

    # Get tntblast results for each assay
    for assay_details in assays:
        run_TNTBLAST(
            assay_details,
            args.g,
            args.outdir,
            args.prefix,
            args.e,
            args.E,
            args.t,
            args.T,
        )
        # tntblast_results = parse_tntblast_output(assay_details, out_path)
        assay_name = assay_details[0]

        # write_amplicons(name + "_" + assay_name, out_path, tntblast_results)

        run_msa(assay_details, args.prefix, args.outdir)

        (msa_results,) = parse_msa_output(
            assay_details, args.prefix, args.outdir, args.g
        )

        write_assay_report(
            assay_details, args.prefix, args.outdir, msa_results, args.cutoff
        )
        write_variants_report(
            assay_details, args.prefix, args.outdir, msa_results, args.cutoff
        )
        print(f"\n {assay_name} Done.\n")


if __name__ == "__main__":
    main()
