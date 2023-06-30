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
import distance
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Blast import NCBIXML
from Bio import SeqIO
__version__ = "0.1.0"
__author__ = "Xiaoli Dong"
pd.set_option("display.float_format", "{:.2f}".format)
iupac_letters = "RYWSMKHBVDN"

iupac_dict = {
    "A": ["A"],
    "G": ["G"],
    "C": ["C"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "W": ["A", "T"],
    "S": ["G", "C"],
    "M": ["A", "C"],
    "K": ["G", "T"],
    "H": ["A", "C", "T"],
    "B": ["C", "G", "T"],
    "V": ["A", "C", "G"],
    "D": ["A", "G", "T"],
    "N": ["A", "C", "G", "T"],
}


def get_parser():

    # Disable default help
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    version = "%(prog)s " + __version__
    parser.add_argument("--version", "-v", action="version", version=version)

    subparsers = parser.add_subparsers(dest="cmd")
    tntblast_parser = subparsers.add_parser(
        "tntblast",
        help="Use tntblast to identify the amplicons from the input sequence file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    tntblast_required_group = tntblast_parser.add_argument_group("required arguments")
    tntblast_optional_group = tntblast_parser.add_argument_group("optional arguments")

    tntblast_required_group.add_argument(
        "--assay",
        "-a",
        type=str,
        required=True,
        help="path to assay details in CSV file",
    )

    # prog give the program basename without path but __file__ gives the full path
    # version = parser.prog + " " + __version__
    # tntblast_optional_group.add_argument(
    #     "--version", "-v", action="version", version=version
    # )

    tntblast_optional_group.add_argument(
        "--outdir",
        "-o",
        type=str,
        default=".",
        help="path to output directory",
    )
    tntblast_optional_group.add_argument(
        "--prefix",
        "-p",
        type=str,
        default="",
        help="output file name prefix",
    )
    tntblast_optional_group.add_argument(
        "--minAbundant",
        type=float,
        default=0,
        help="minimum prevalence (%%) of primer site variants reported in reports, min > 0, max < 100",
    )
    tntblast_optional_group.add_argument(
        "-e",
        type=int,
        default=50,
        help="The minimum allowed temperature (in °C) for a primer oligo to bind to a target sequence",
    )
    tntblast_optional_group.add_argument(
        "-E",
        type=int,
        default=50,
        help="The minimum allowed temperature (in °C) for a probe oligo to bind to a target sequence",
    )
    tntblast_optional_group.add_argument(
        "-t",
        type=float,
        default=0.90,
        help="The Molar concentration of primer oligos, used to compute double stranded DNA melting temperatures",
    )
    tntblast_optional_group.add_argument(
        "-T",
        type=float,
        default=0.20,
        help="The Molar concentration of probe oligos, used to compute double stranded DNA melting temperatures",
    )
    tntblast_optional_group.add_argument(
        "--mask_output",
        default=False,
        action="store_true",
        help="Turn on output mask options to convert the identical bases as dot ",
    )
    tntblast_optional_group.add_argument(
        "--enable_iupac",
        default=False,
        action="store_true",
        help="Turn on iupac option to enable iupac degeneration during the sequence mask ",
    )
    tntblast_optional_group.add_argument(
        "--exclude_ns_from_amplicon",
        default=False,
        action="store_true",
        help="Turn on this option will exclude the amplicons containing Ns",
    )
    tntblast_parser.set_defaults(func=run_tntblast_analysis)

    blastn_parser = subparsers.add_parser(
        "blastn",
        help="Use NCBI blastn to identify the amplicons from the input sequence file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    blastn_required_group = blastn_parser.add_argument_group("required arguments")
    blastn_optional_group = blastn_parser.add_argument_group("optional arguments")

    blastn_required_group.add_argument(
        "--assay",
        "-a",
        type=str,
        required=True,
        help="path to assay details in CSV file",
    )
    blastn_required_group.add_argument(
        "--template", type=str, required=True, help="fasta format PCR template file"
    )

    blastn_optional_group.add_argument(
        "--outdir",
        "-o",
        type=str,
        default=".",
        help="path to output directory",
    )
    blastn_optional_group.add_argument(
        "--prefix",
        type=str,
        default="",
        help="output file name prefix",
    )

    blastn_optional_group.add_argument(
        "--minAbundant",
        type=float,
        default=0,
        help="minimum prevalence (%%) of primer site variants reported in reports, min > 0, max < 100",
    )
    blastn_optional_group.add_argument(
        "--word_size",
        type=int,
        default=7,
        help="blastn search word size for wordfinder algorithm (length of best perfect match)",
    )

    blastn_optional_group.add_argument(
        "--evalue",
        type=float,
        default=10,
        help="blastn search expectation value (E) threshold for saving hits",
    )
    blastn_optional_group.add_argument(
        "--perc_identity",
        type=float,
        default=80,
        help="blastn search Percent identity between 0..100",
    )
    blastn_optional_group.add_argument(
        "--qcov_hsp_perc",
        type=float,
        default=95,
        help="blastn search Percent query coverage per hsp between 0..100",
    )
    
    blastn_optional_group.add_argument(
        "--max_target_seqs",
        "-m",
        type=int,
        default=1000000,
        help="blastn maximum number of aligned sequences to keep for each search",
    )
    # Switch
    blastn_optional_group.add_argument(
        "--subject_besthit",
        default=False,
        action="store_true",
        help="Turn on best hit per subject sequence for blastn search",
    )
    blastn_optional_group.add_argument(
        "--exclude_ns_from_amplicon",
        default=False,
        action="store_true",
        help="Turn on this option will exclude the amplicons containing Ns",
    )
    

    blastn_optional_group.add_argument(
        "--mask_output",
        default=False,
        action="store_true",
        help="Turn on output mask options to convert the identical bases as dot ",
    )
    blastn_optional_group.add_argument(
        "--enable_iupac",
        default=False,
        action="store_true",
        help="Turn on iupac option to enable iupac degeneration during the sequence mask ",
    )
    blastn_parser.set_defaults(func=run_blastn_analysis)

    return parser


def set_pathes(outdir, prefix, assay_name, search_tool_name):

    if prefix:
        path_to_assay_seq = os.path.join(
            outdir + "/tmp", prefix + "_" + assay_name + "_msa_assay.fasta"
        )
        path_to_amplicon_before_correction = os.path.join(
            outdir + "/tmp", 
            prefix + "_" + assay_name + "_" + search_tool_name + "_amplicon.no_correction.fasta",
        )
        path_to_amplicon = os.path.join(
            outdir,
            prefix + "_" + assay_name + "_" + search_tool_name + "_amplicon.fasta",
        )
        path_to_amplicon_exclude = os.path.join(
            outdir + "/tmp", 
            prefix + "_" + assay_name + "_" + search_tool_name + "_amplicon.exclude.fasta",
        )
        
        path_to_mafft_output = os.path.join(
            outdir,
            prefix + "_" + assay_name + "_" + search_tool_name + "_amplicon_mafft.fasta",
        )

        path_to_pcr_report = os.path.join(
            outdir,
            prefix + "_" + assay_name + "_" + search_tool_name + "_pcr.tsv",
        )

        path_to_assay_report = os.path.join(
            outdir,
            prefix + "_" + assay_name + "_" + search_tool_name + "_assay_report.tsv",
        )        
    else:
        path_to_assay_seq = os.path.join(
            outdir + "/tmp", assay_name + "_msa_assay.fasta"
        )
        path_to_amplicon_before_correction = os.path.join(
            outdir + "/tmp", 
           assay_name + "_" + search_tool_name + "_amplicon.no_correction.fasta",
        )
        path_to_amplicon = os.path.join(
            outdir,
           assay_name + "_" + search_tool_name + "_amplicon.fasta",
        )
        path_to_amplicon_exclude = os.path.join(
            outdir + "/tmp", 
           assay_name + "_" + search_tool_name + "_amplicon.exclude.fasta",
        )
        
        path_to_mafft_output = os.path.join(
            outdir,
           assay_name + "_" + search_tool_name + "_amplicon_mafft.fasta",
        )

        path_to_pcr_report = os.path.join(
            outdir,
           assay_name + "_" + search_tool_name + "_pcr.tsv",
        )

        path_to_assay_report = os.path.join(
            outdir,
           assay_name + "_" + search_tool_name + "_assay_report.tsv",
        )
        
    return (
        path_to_assay_seq,
        path_to_amplicon_before_correction,
        path_to_amplicon,
        path_to_amplicon_exclude,
        path_to_mafft_output,
        path_to_pcr_report,
        path_to_assay_report,
        
    )


def run_blastn_analysis(args):

    # check output dir exists otherwise create
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    if args.minAbundant < 0 or args.minAbundant > 100:
        print(
            f"\nERROR: minimum prevalence (%) of primer site variants should in the range [0,100]\n"
        )
        exit(1)

    # Check input genomes file
    # check_genomes_file(args.g)

    # Parse assay file to get oligo names and seqs
    assays = read_assay_file(args.assay)

    search_tool_name = "blastn"
    for assay_details in assays:
        assay_name = assay_details[0]
        fwd_primer_name = assay_details[2]
        fwd_primer_seq = assay_details[3]
        rev_primer_name = assay_details[4]
        rev_primer_seq = assay_details[5]
        probe_name = assay_details[6]
        probe_seq = assay_details[7]

        ref_seq_file_path = assay_details[1]
        assay_outdir = os.path.join(args.outdir, assay_name)
        Path(assay_outdir + "/tmp").mkdir(parents=True, exist_ok=True)
        (
            path_to_assay_seq,
            path_to_amplicon_before_correction,
            path_to_amplicon,
            path_to_amplicon_exclude,
            path_to_mafft_output,
            path_to_pcr_report,
            path_to_assay_report,
            
        ) = set_pathes(assay_outdir, args.prefix, assay_name, search_tool_name)

        # # the orientiation of the sequences is on the plus strand

        f = open(path_to_assay_seq, "w", encoding="utf-8")
        f.write(f">{fwd_primer_name}\n{fwd_primer_seq}\n")
        f.write(
            f">{rev_primer_name}\n{str(Seq(rev_primer_seq).reverse_complement())}\n"
        )
        if probe_name != "" and probe_seq != "":
            f.write(f">{probe_name}\n{probe_seq}\n")
        f.close()

        # path_to_blastdb = run_makeblastdb(args.g)
        path_to_blastdb = run_makeblastdb(ref_seq_file_path, assay_name)

        path_to_blastn_xml_output = run_blastn(
            args.template,
            path_to_blastdb,
            assay_details,
            args.prefix,
            args.outdir,
            args.max_target_seqs,
            args.word_size,
            args.evalue,
            args.subject_besthit,
            args.perc_identity,
            args.qcov_hsp_perc
        )
        #path_to_amplicon,amplicon_coordinates  = parse_blastn_output(path_to_blastn_xml_output, path_to_amplicon, args.qcov)
        path_to_amplicon = parse_blastn_output(path_to_blastn_xml_output, args.exclude_ns_from_amplicon, path_to_amplicon_before_correction, path_to_amplicon, path_to_amplicon_exclude, ref_seq_file_path)
        
        #refine_amplicon(ref_seq_file_path, path_to_amplicon, amplicon_coordinates)
        #print(amplicon_coordinates)
        if os.stat(path_to_amplicon).st_size == 0:
            print(
                f"there is no identified target amplicons for {assay_name} using blastn and proceed to the next assay ... "
            )
            continue

        run_msa(path_to_assay_seq, path_to_amplicon, path_to_mafft_output)
        (msa_results, my_tech_dict) = parse_msa_output(
            assay_details,
            ref_seq_file_path,
            path_to_amplicon,
            path_to_mafft_output,
            path_to_pcr_report,
        )

        write_assay_report(
            assay_details,
            msa_results,
            args.minAbundant,
            args.mask_output,
            args.enable_iupac,
            path_to_assay_report,
            my_tech_dict,
        )
        


def run_tntblast_analysis(args):

    # check output dir exists otherwise create
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    if args.minAbundant < 0 or args.minAbundant > 100:
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

    assays = read_assay_file(args.assay)

    # Create empty dataframe for tabulated tntblast results

    # Get tntblast results for each assay
    search_tool_name = "tntblast"
    for assay_details in assays:
        assay_name = assay_details[0]
        fwd_primer_name = assay_details[2]
        fwd_primer_seq = assay_details[3]
        rev_primer_name = assay_details[4]
        rev_primer_seq = assay_details[5]
        probe_name = assay_details[6]
        probe_seq = assay_details[7]
        ref_seq_file_path = assay_details[1]
        assay_outdir = os.path.join(args.outdir, assay_name)
        Path(assay_outdir + "/tmp").mkdir(parents=True, exist_ok=True)
        #Path(assay_outdir).mkdir(parents=True, exist_ok=True)
        (
            path_to_assay_seq,
            path_to_amplicon_before_correction,
            path_to_amplicon,
            path_to_amplicon_exclude,
            path_to_mafft_output,
            path_to_pcr_report,
            path_to_assay_report,
        ) = set_pathes(assay_outdir, args.prefix, assay_name, search_tool_name)

        f = open(path_to_assay_seq, "w", encoding="utf-8")
        f.write(f">{fwd_primer_name}\n{fwd_primer_seq}\n")
        f.write(
            f">{rev_primer_name}\n{str(Seq(rev_primer_seq).reverse_complement())}\n"
        )
        if probe_name != "" and probe_seq != "":
            f.write(f">{probe_name}\n{probe_seq}\n")
        f.close()

        run_TNTBLAST(
            args.exclude_ns_from_amplicon,
            assay_details,
            ref_seq_file_path,
            assay_outdir,
            args.e,
            args.E,
            args.t,
            args.T,
            path_to_amplicon_before_correction,
            path_to_amplicon,
            path_to_amplicon_exclude
        )
        if os.stat(path_to_amplicon).st_size == 0:
            print(
                f"there is no identified target amplicons for {assay_name} and proceed to the next assay ... "
            )
            continue

        run_msa(path_to_assay_seq, path_to_amplicon, path_to_mafft_output)

        (msa_results, my_tech_dict) = parse_msa_output(
            assay_details,
            ref_seq_file_path,
            path_to_amplicon,
            path_to_mafft_output,
            path_to_pcr_report,
        )

        write_assay_report(
            assay_details,
            msa_results,
            args.minAbundant,
            args.mask_output,
            args.enable_iupac,
            path_to_assay_report,
            my_tech_dict,
        )
        

def check_genomes_file(path_to_file, assay_name):
    """Checks that file exists, is not empty, and headers are unique and hashable."""
    if os.path.exists(path_to_file) == False or os.path.isfile(path_to_file) == False:
        print(
            f"\nERROR: Genomes FASTA file {path_to_file} for {assay_name} does not exist.\n"
        )
        exit(1)
    if os.path.getsize(path_to_file) == 0:
        print(
            f"\nERROR: Genomes FASTA file {path_to_file}  for {assay_name}  is empty.\n"
        )
        exit(1)
    with open(path_to_file, "r", encoding="utf-8") as input_file:
        headers = []
        for line in input_file:
            if line[0] == ">":
                headers.append(line.strip().lstrip(">"))
        if len(headers) != len(set(headers)):
            print(
                "\nERROR: {assay_name} FASTA headers in input genomes file must be unique.\n"
            )
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

    with open(path_to_file, "r", encoding="utf-8") as input_file:
        for line in input_file:
            if line.startswith("#"):
                continue
            line = line.rstrip().split(",")
            # Check if line has the accepted number of comma-separated values
            if len(line) not in [6, 8]:
                print("\nERROR: Line in assay file not properly formatted:")
                print(",".join(line))
                print("\nAssay file lines must be a comma-separated list:")
                print(
                    "assay_name,reference_seq_file_path,fwd_primer_name,fwd_primer_seq,rev_primer_name,rev_primer_seq,probe_name,probe_seq"
                )
                exit(1)
            else:

                assay_name = line[0]
                ref_seq_file_path = line[1]
                fwd_primer_name = line[2]
                fwd_primer_seq = line[3]
                rev_primer_name = line[4]

                rev_primer_seq = line[5]
                if len(line) == 8:
                    probe_name = line[6]
                    probe_seq = line[7]
                else:
                    probe_name = ""
                    probe_seq = ""
                check_genomes_file(ref_seq_file_path, assay_name)
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

                        print()
                        exit(1)
                # Create assay info tuple and append to assays list
                assay = (
                    assay_name,
                    ref_seq_file_path,
                    fwd_primer_name,
                    fwd_primer_seq,
                    rev_primer_name,
                    rev_primer_seq,
                    probe_name,
                    probe_seq,
                )

                assays.append(assay)
    # Check that names used for assays and oligos are unique
    for name_index in [0, 2, 4, 6]:
        names = [assay[name_index] for assay in assays if assay[name_index] != ""]
        if len(names) != len(set(names)):
            print("\nERROR: Assay and oligo names must be unique!\n")
            exit()
    return tuple(assays)


def run_TNTBLAST(
    exclude_ns_from_amplicon,
    assay_details,
    path_to_genomes,
    outdir,
    melting_temp_primer,
    melting_temp_probe,
    primer_molarity,
    probe_molarity,
    path_to_amplicon_before_correction,
    path_to_amplicon,
    path_to_amplicon_exclude
):
    """Run TNTBLAST with provided assay details."""

    # Create TNTBLAST assay file from assay details (TSV file of assay name and oligo seqs)
    assay_name = assay_details[0]
    fwd_primer_seq = assay_details[3]
    rev_primer_seq = assay_details[5]
    probe_seq = assay_details[7]
    print(f"Running TNTBLAST for {assay_name} against {path_to_genomes}...")
    print(f"Fwd primer seq: {fwd_primer_seq}\nRev primer seq: {rev_primer_seq}")
    if probe_seq != "":
        print(f"Probe seq: {probe_seq}")
    assay = [assay_name, fwd_primer_seq, rev_primer_seq]
    assay += [probe_seq] if probe_seq != "" else []
    assay = "\t".join(assay)
    path_to_tntblast_assay = os.path.join(outdir + "/tmp", assay_name + "_tntblast_assay.tsv")
    with open(path_to_tntblast_assay, "w", encoding="utf-8") as output_file:
        output_file.write(assay + "\n")
    # Create terminal command for TNTBLAST and run
    primer_molarity = primer_molarity / 1000000
    probe_molarity = probe_molarity / 1000000
    terminal_command = (
        f"tntblast -i {path_to_tntblast_assay} -d {path_to_genomes} -o {path_to_amplicon_before_correction}"
        f" -e {melting_temp_primer} -E {melting_temp_probe} -t {primer_molarity:.2e} -T {probe_molarity:.2e}"
        f" --best-match -m 1 -v F --single-primer-pcr F"
    )
    print(terminal_command)
    completed_process = subprocess.run(
        terminal_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True,
        check=True,
    )
    if completed_process.returncode != 0:
        print(
            f"\nERROR: TNTBLAST terminated with errors.\nTNTBLAST error code: {completed_process.returncode}\n"
        )
        exit(1)
    # GARBAGE COLLECTION (remove TNTBLAST assay TSV file)
    # os.remove(path_to_tntblast_assay)
    print()
    # return path_to_tntblast_txt_output
    with open(path_to_amplicon,"w") as f, open(path_to_amplicon_exclude, "w") as f0:

        for amplicon_record in SeqIO.parse(path_to_amplicon_before_correction, "fasta"):
            if exclude_ns_from_amplicon:
                if "N".lower() not in str(amplicon_record.seq).lower():
                    SeqIO.write(amplicon_record, f, "fasta-2line")
                else:
                    SeqIO.write(amplicon_record, f0, "fasta-2line")
            else:
                SeqIO.write(amplicon_record, f, "fasta-2line")

    f.close()
    f0.close()
    return path_to_amplicon



def run_makeblastdb(path_to_genomes, assay_name):

    """make blastn search database"""
    print(f"Making a blastn search database using {path_to_genomes}...")

    genome_name = os.path.basename(path_to_genomes)

    path_to_blastdb = f"tmp/{assay_name}/{genome_name}"
    terminal_command = (
        f"makeblastdb -in {path_to_genomes} -out {path_to_blastdb} -dbtype nucl"
    )

    print(terminal_command)
    completed_process = subprocess.run(
        terminal_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True,
        check=True,
    )

    if completed_process.returncode != 0:
        print(
            f"\nERROR: makeblastdb terminated with errors.\nmakeblastdb error code: {completed_process.returncode}\n"
        )
        exit(1)
    # GARBAGE COLLECTION (remove TNTBLAST assay TSV file)
    # os.remove(path_to_tntblast_assay)
    print()
    return path_to_blastdb


def run_blastn(
    path_to_query,
    path_to_blastdb,
    assay_details,
    prefix,
    outdir,
    max_target_seqs,
    word_size,
    evalue,
    subject_besthit,
    perc_identity,
    qcov_hsp_perc
):

    """make blastn search database"""
    print(f"Running blastn search against {path_to_blastdb}...")

    # genome_name = os.path.basename(path_to_genomes)
    assay_name = assay_details[0]
    path_to_blastn_xml_output = os.path.join(
        outdir + "/" + assay_name + "/tmp", prefix + "_" + assay_name + "_blastn.xml"
    )
    # path_to_blastdb = f"tmp/{genome_name}"
    if subject_besthit:
        terminal_command = f"blastn -query {path_to_query} -db {path_to_blastdb} -word_size {word_size} -out {path_to_blastn_xml_output}  -outfmt 5  -max_target_seqs {max_target_seqs}  -dust no -task blastn -evalue {evalue} -perc_identity {perc_identity} -qcov_hsp_perc {qcov_hsp_perc} -subject_besthit"
    else:
        terminal_command = f"blastn -query {path_to_query} -db {path_to_blastdb} -word_size {word_size} -out {path_to_blastn_xml_output}  -outfmt 5  -max_target_seqs {max_target_seqs}  -dust no -task blastn -evalue {evalue} -perc_identity {perc_identity} -qcov_hsp_perc {qcov_hsp_perc}"
    print(terminal_command)

    completed_process = subprocess.run(
        terminal_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True,
        check=True,
    )

    if completed_process.returncode != 0:
        print(
            f"\nERROR: blastn terminated with errors.\nblastn error code: {completed_process.returncode}\n"
        )
        exit(1)
    # GARBAGE COLLECTION (remove TNTBLAST assay TSV file)
    # os.remove(path_to_tntblast_assay)
    print()
    return path_to_blastn_xml_output


# def parse_blastn_output(assay_details, path_to_blastn_xml, outdir, prefix):
def parse_blastn_output(path_to_blastn_xml, exclude_ns_from_amplicon, path_to_amplicon_before_correction, path_to_amplicon, path_to_amplicon_exclude, ref_seq_file_path):

    # https://uoftcoders.github.io/studyGroup/lessons/python/biopython/lesson/
    # path_to_amplicon = os.path.join(
    #     outdir, prefix + "_" + assay_name + "_blastn_amplicon.fasta"
    # )
    f = open(path_to_amplicon_before_correction, "w", encoding="utf-8")

    result_handle = open(path_to_blastn_xml, "r")
    blast_records = NCBIXML.parse(result_handle)

    # get rid of the duplicate ids when there are multiple input templates
    seen = []
    amplicon_coordinates = {}

    for blast_record in blast_records:
        qlen = blast_record.query_letters
        
        for alignment in blast_record.alignments:
           
            for hsp in alignment.hsps:
                #if hsp.align_length / qlen >= qcov:
                # print("****Alignment****")
                hit_def = alignment.hit_def
                id = hit_def.split(" ", 1)[0]
                
                if hsp.strand[1]:
                    coord = {
                        "sbjct_start":hsp.sbjct_start, 
                        "sbjct_end":hsp.sbjct_end,
                        "query_start": hsp.query_start,
                        "query_end": hsp.query_end,
                        "query_len": qlen,
                        "strand": hsp.strand[1],
                        "end2end": True if qlen == hsp.query_end else False
                        }
                else:
                    coord = {
                    "sbjct_start":hsp.sbjct_end, 
                    "sbjct_end":hsp.sbjct_start,
                    "query_start": hsp.query_start,
                    "query_end": hsp.query_end,
                    "query_len": qlen,
                    "strand": hsp.strand[1],
                    "end2end": True if qlen == hsp.query_start else False
                    }

                if id not in seen:
                    amplicon_coordinates[id] = coord
                    seen.append(id)
                    #f.write(f">{id}\n{hsp.sbjct}\n")
                    f.write(f">{hit_def}\n{hsp.sbjct}\n")
    f.close()

    refine_amplicon(exclude_ns_from_amplicon, ref_seq_file_path, path_to_amplicon_before_correction, path_to_amplicon, path_to_amplicon_exclude, amplicon_coordinates)
    #return path_to_amplicon, amplicon_coordinates
    return path_to_amplicon

def refine_amplicon(exclude_ns_from_amplicon, path_to_genomes, path_to_amplicon_before_correction, path_to_amplicon, path_to_amplicon_exclude, amplicon_coordinates):

    SeqDict = SeqIO.to_dict(SeqIO.parse(path_to_genomes, "fasta"))
    print(path_to_amplicon)
    with open(path_to_amplicon,"w") as f, open(path_to_amplicon_exclude, "w") as f0:
        for amplicon_record in SeqIO.parse(path_to_amplicon_before_correction, "fasta"):
                #f.write(str(amplicon_record.id) + "\n")
                # f.write(str(amplicon_record.seq[:10]) + "\n")  #first 10 base positions
                # f.write(str(amplicon_record.seq[-10:]) + "\n") #last 10 base positions
                id = amplicon_record.id
                ref_record = SeqDict[id]
                s_start = amplicon_coordinates[id]["sbjct_start"]
                s_end = amplicon_coordinates[id]["sbjct_end"]
                query_start = amplicon_coordinates[id]["query_start"]
                query_end = amplicon_coordinates[id]["query_end"]
                
                qlen = amplicon_coordinates[id]["query_len"]
                s_strand = amplicon_coordinates[id]["strand"]
                t_start = s_start
                t_end = s_end
                #print(f"id={id}, start={s_start}, end={s_end}, strand={s_strand}")
                if amplicon_coordinates[amplicon_record.id]["end2end"]:
                        if exclude_ns_from_amplicon:
                            if "N".lower() not in str(amplicon_record.seq).lower(): 
                                SeqIO.write(amplicon_record, f, "fasta-2line")
                            else:
                                SeqIO.write(amplicon_record, f0, "fasta-2line")
                        else:
                            SeqIO.write(amplicon_record, f, "fasta-2line")
                        
                else:
                    if s_strand == "Plus":
                        if query_start > 1:
                            t_start = t_start - query_start
                        if query_end < qlen:
                            t_end = t_end + qlen - query_end
                        # blast coordinate with 1 but python index starting with 0
                        cut_record = ref_record[t_start-1: t_end]
                        cut_record.id = ref_record.id
                        cut_record.name = ""
                        cut_record.description = ref_record.description
                        if exclude_ns_from_amplicon:
                            if "N".lower() not in str(cut_record.seq).lower(): 
                                SeqIO.write(cut_record, f, "fasta-2line")
                            else:
                                SeqIO.write(cut_record, f0, "fasta-2line")
                        else:
                            SeqIO.write(cut_record, f, "fasta-2line")
                       
                    else:
                        #print("minus.........................")
                        if query_start > 1:
                            # blast coordinate with 1 but python index starting with 0
                            t_start = t_start + query_start
                        if query_end < qlen:
                            t_end = t_end - qlen + query_end
                        
                        #cut_record = ref_record[t_end-1: t_start].reverse_complement()
                        #cut_record = ref_record[t_end-1: t_start].reverse_complement(id=f"{ref_record.id}_{t_start}_{t_end}_{s_strand}",  description="")
                        cut_record = ref_record[t_end-1: t_start].reverse_complement(id=f"{ref_record.id}_rc", name=True, description=True)
                        if exclude_ns_from_amplicon:
                            if "N".lower() not in str(amplicon_record.seq).lower(): 
                                SeqIO.write(cut_record, f, "fasta-2line")
                            else:
                                SeqIO.write(cut_record, f0, "fasta-2line")
                        else:
                            SeqIO.write(cut_record, f, "fasta-2line")

    f.close()
    f0.close()             

def count_targets(path_to_genomes):
    """Counts the number of target sequences in the provided genomes FASTA file.
    Returns an int of the count."""
    with open(path_to_genomes, "r", encoding="utf-8") as input_file:
        seq_counter = 0
        for line in input_file:
            if line[0] == ">":
                seq_counter += 1
    return int(seq_counter)


def get_targets(path_to_genomes):
    """Reads all sequence in the provided genomes FASTA file and returns them in
    a dict where keys are FASTA headers and values are the nucleotide sequences."""
    with open(path_to_genomes, "r", encoding="utf-8") as input_file:
        seqs = {}
        for line in input_file:
            if line[0] == ">":
                header = line.strip().lstrip(">")
                seqs[header] = ""
            else:
                seqs[header] += line.strip()
    return seqs


# def run_msa(assay_details, prefix, outdir, search_tool_name):
def run_msa(path_to_assay_seq, path_to_amplicon, path_to_mafft_output):

    """Run mafft with provided assay details."""

    print(
        f"Running mafft to produce amplicon multiple sequence alignment for {path_to_amplicon}..."
    )

    # terminal_command = f"cat {path_to_assay_seq} {path_to_amplicon} | mafft --auto --thread -1 - > {path_to_mafft_output}"
    terminal_command = f"mafft --auto --thread -1 {path_to_amplicon} | mafft --auto --thread -1 --addfragments {path_to_assay_seq} - > {path_to_mafft_output}"
    print(terminal_command)
    completed_process = subprocess.run(
        terminal_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True,
        check=True,
    )
    if completed_process.returncode != 0:
        print(
            f"\nERROR: mafft terminated with errors.\nmafft error code: {completed_process.returncode}\n"
        )
        exit(1)
    # GARBAGE COLLECTION (remove TNTBLAST assay TSV file)
    # os.remove(path_to_tntblast_assay)
    print()


# def parse_msa_output(assay_details, prefix, outdir, path_to_genomes, search_tool_name):
def parse_msa_output(
    assay_details, path_to_genomes, path_to_amplicon, path_to_mafft, path_to_pcr_report
):
    assay_name = assay_details[0]
    fwd_primer_name = assay_details[2]
    rev_primer_name = assay_details[4]
    probe_name = assay_details[6]

    ################################

    total_input_seqs = count_targets(path_to_genomes)
    amplicon_count = count_targets(path_to_amplicon)
    #################################

    align = AlignIO.read(path_to_mafft, "fasta")

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
    tech_dict = {}
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
            #for record in align[:, 0:fwd_end]:
            for record in align[:, fwd_start:fwd_end]:
                # print(str(fwd_record.seq))
                # if "n" in str(record.seq):
                #     continue
                if record.id == fwd_primer_name:
                    tech_dict["fwd_primer_align"] = str(record.seq)
                    continue
                elif (record.id == rev_primer_name) or (record.id == probe_name):
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
            #for record in align[:, rev_start : len(rev_primer_alignment)]:
            for record in align[:, rev_start : rev_end]:
                # for record amplicon_recordalign[:, rev_start:rev_end]:
                # print(str(fwd_record.seq))
                # if "n" in str(record.seq):
                #     continue
                if record.id == rev_primer_name:
                    tech_dict["rev_primer_align"] = str(record.seq)
                    continue
                elif (record.id == fwd_primer_name) or (record.id == probe_name):
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
                # if "n" in str(record.seq):
                #      continue
                
                if record.id == probe_name:
                    tech_dict["probe_align"] = str(record.seq)
                    continue
                elif (record.id == fwd_primer_name) or (record.id == rev_primer_name):
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

        if assay_details[6] == "" and fwd_coord_found and rev_coord_found:
            # msa_results = pd.merge(fwd_df_mask, rev_df_mask, on="target")
            msa_results = pd.merge(fwd_df, rev_df, on="target")
            msa_results["total_errors"] = (
                msa_results["fwd_errors"] + msa_results["rev_errors"]
            )
            msa_results = msa_results[
                (msa_results.target != fwd_primer_name)
                & (msa_results.target != rev_primer_name)
            ]

            """ msa_results.to_csv(
                path_to_pcr_report,
                sep="\t",
                index=False,
                # index_label="variant_name",
            )
            """
            # print(msa_results)
            break
        elif (
            assay_details[6] != ""
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

            """ msa_results.to_csv(
                path_to_pcr_report,
                sep="\t",
                index=False,
                # index_label="variant_name",
            ) """

            break
    print(tech_dict)
    msa_results["assay_name"] = assay_name
    msa_results["total_input_seqs"] = total_input_seqs
    msa_results["amplicon_detected"] = amplicon_count
    # print(msa_results)
    return (msa_results, tech_dict)


def write_assay_report(
    assay_details,
    msa_results,
    threshold,
    mask_output,
    enable_iupac,
    path_to_assay_report,
    my_tech_dict,
):
    total_input_seqs = msa_results["total_input_seqs"].iloc[0]
    amplicon_count = msa_results["amplicon_detected"].iloc[0]
    # taqman assay with probe
    if assay_details[6] != "":
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
            msa_results.groupby(all_cols).size().reset_index(name="amplicon_detected")
        )

        all_summary["total_count_seqdb"] = total_input_seqs
        all_summary["total_amplicon"] = amplicon_count
        all_summary["total_amplion_of_seqdb_pct"] = round(
            amplicon_count * 100 / total_input_seqs, 2
        )

        all_summary["amplicon_detected_of_total_amplicon_pct"] = round(
            all_summary["amplicon_detected"] * 100 / amplicon_count, 2
        )

        all_summary.astype(
            {
                "total_errors": "int",
                "amplicon_detected": "int",
                "total_amplicon": "int",
                "total_count_seqdb": "int",
            }
        )
        all_summary = all_summary[
            all_summary["amplicon_detected_of_total_amplicon_pct"] >= float(threshold)
        ]
        # all_summary.sort_values(
        #     ["amplicon_detected"], axis=0, ascending=[False], inplace=True
        # )
        # print(all_summary)
        new_cols = [
            "fwd_site",
            "fwd_errors",
            "probe_site",
            "probe_errors",
            "rev_site",
            "rev_errors",
            "total_errors",
            "total_count_seqdb",
            "total_amplicon",
            "total_amplion_of_seqdb_pct",
            "amplicon_detected",
            "amplicon_detected_of_total_amplicon_pct",
        ]
        all_summary = all_summary.reindex(columns=new_cols)
        # print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        # print(all_summary)
        ################################### start mask ####################
        if mask_output:
            fwd_list = list(my_tech_dict["fwd_primer_align"])
            index = 0
            for seqstr in all_summary["fwd_site"]:

                # create a lit using the sequence string
                seqstr_list = list(seqstr)
                seqstr_list_masked = []
                errors = 0
                for i, letter in enumerate(seqstr_list):
                    if enable_iupac:
                        if letter.upper() in iupac_dict[fwd_list[i].upper()]:
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)
                            errors += 1
                    else:
                        if letter == fwd_list[i]:
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)

                # update the seq str with the masked str, dot represent the same as the input tech seq
                all_summary.at[index, "fwd_site"] = "".join(seqstr_list_masked)
                if enable_iupac:
                    all_summary.at[index, "fwd_errors"] = errors
                index += 1

                # print(seqstr_list_masked)
                # print(all_summary)

            rev_list = list(my_tech_dict["rev_primer_align"])
            index = 0

            for seqstr in all_summary["rev_site"]:
                errors = 0
                # create a lit using the sequence string
                seqstr_list = list(seqstr)
                seqstr_list_masked = []
                for i, letter in enumerate(seqstr_list):
                    if enable_iupac:
                        if letter.upper() in iupac_dict[rev_list[i].upper()]:
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)
                            errors += 1
                    else:
                        if letter == rev_list[i]:
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)

                # update the seq str with the masked str, dot represent the same as the input tech seq
                all_summary.at[index, "rev_site"] = "".join(seqstr_list_masked)
                if enable_iupac:
                    all_summary.at[index, "rev_errors"] = errors
                index += 1

            probe_list = list(my_tech_dict["probe_align"])
            index = 0
            for seqstr in all_summary["probe_site"]:
                # print("seqstr=" + seqstr + ", probe=" + my_tech_dict["probe_align"])
                errors = 0
                # create a lit using the sequence string
                seqstr_list = list(seqstr)
                seqstr_list_masked = []
                for i, letter in enumerate(seqstr_list):
                    if enable_iupac:
                        if letter.upper() in iupac_dict[probe_list[i].upper()]:
                            # print("letter=" + letter)
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)
                            errors += 1
                    else:
                        if letter == probe_list[i]:
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)

                # update the seq str with the masked str, dot represent the same as the input tech seq
                # print("".join(seqstr_list_masked))
                all_summary.at[index, "probe_site"] = "".join(seqstr_list_masked)
                if enable_iupac:
                    all_summary.at[index, "probe_errors"] = errors
                index += 1
        ################################# end of mask ###########################

        # if enable_iupac:
        # update total_errors
        all_summary["total_errors"] = (
            all_summary["fwd_errors"]
            + all_summary["rev_errors"]
            + all_summary["probe_errors"]
        )
        # print(all_summary)
        # when iupac enable, after degenerated, some of the different site will become the same
        final_all_summary = all_summary.groupby(
            [
                "fwd_site",
                "fwd_errors",
                "probe_site",
                "probe_errors",
                "rev_site",
                "rev_errors",
                "total_errors",
                "total_count_seqdb",
                "total_amplicon",
                "total_amplion_of_seqdb_pct",
            ]
        ).sum()

        final_all_summary.sort_values(
            ["amplicon_detected"], axis=0, ascending=[False], inplace=True
        )
        final_all_summary = final_all_summary.reset_index()

        # print(final_all_summary)
        final_all_summary["amplicon_detected_of_total_amplicon_pct"] = round(
            final_all_summary["amplicon_detected_of_total_amplicon_pct"], 2
        )
        final_all_summary["total_amplion_of_seqdb_pct"] = round(
            final_all_summary["total_amplion_of_seqdb_pct"], 2
        )

        final_all_summary.to_csv(
            path_to_assay_report,
            sep="\t",
            index=True,
            # index_label="variant_id",
        )
    else:
        # print(msa_results)
        all_cols = ["fwd_site", "fwd_errors", "rev_site", "rev_errors", "total_errors"]
        all_summary = (
            msa_results.groupby(all_cols).size().reset_index(name="amplicon_detected")
        )

        all_summary["total_count_seqdb"] = total_input_seqs
        all_summary["total_amplicon"] = amplicon_count
        all_summary["total_amplion_of_seqdb_pct"] = round(
            amplicon_count * 100 / total_input_seqs, 2
        )
        all_summary["amplicon_detected_of_total_amplicon_pct"] = round(
            all_summary["amplicon_detected"] * 100 / amplicon_count, 2
        )

        all_summary.astype(
            {
                "total_errors": "int",
                "amplicon_detected": "int",
                "total_amplicon": "int",
                "total_count_seqdb": "int",
            }
        )
        all_summary = all_summary[
            all_summary["amplicon_detected_of_total_amplicon_pct"] >= float(threshold)
        ]

        new_cols = [
            "fwd_site",
            "fwd_errors",
            "rev_site",
            "rev_errors",
            "total_errors",
            "total_count_seqdb",
            "total_amplicon",
            "total_amplion_of_seqdb_pct",
            "amplicon_detected",
            "amplicon_detected_of_total_amplicon_pct",
        ]
        all_summary = all_summary.reindex(columns=new_cols)

        ################################### start mask ####################
        if mask_output:
            fwd_list = list(my_tech_dict["fwd_primer_align"])
            index = 0
            for seqstr in all_summary["fwd_site"]:
                errors = 0
                # create a lit using the sequence string
                seqstr_list = list(seqstr)
                seqstr_list_masked = []
                for i, letter in enumerate(seqstr_list):

                    if enable_iupac:
                        if letter.upper() in iupac_dict[fwd_list[i].upper()]:
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)
                            errors += 1
                    else:
                        if letter == fwd_list[i]:
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)

                # update the seq str with the masked str, dot represent the same as the input tech seq
                all_summary.at[index, "fwd_site"] = "".join(seqstr_list_masked)

                index += 1

                # print(seqstr_list_masked)
                # print(all_summary)

            rev_list = list(my_tech_dict["rev_primer_align"])
            index = 0
            for seqstr in all_summary["rev_site"]:
                errors = 0
                # create a lit using the sequence string
                seqstr_list = list(seqstr)
                seqstr_list_masked = []

                for i, letter in enumerate(seqstr_list):
                    if enable_iupac:
                        if letter.upper() in iupac_dict[rev_list[i].upper()]:
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)
                            errors += 1
                    else:
                        if letter == rev_list[i]:
                            seqstr_list_masked.append(".")
                        else:
                            seqstr_list_masked.append(letter)

                # update the seq str with the masked str, dot represent the same as the input tech seq
                all_summary.at[index, "rev_site"] = "".join(seqstr_list_masked)
                index += 1

                # print(seqstr_list_masked)
                # print(all_summary)
        ##############################
        # print(all_summary)

        if enable_iupac:
            # update total_errors
            all_summary["total_errors"] = (
                all_summary["fwd_errors"]
                + all_summary["rev_errors"]
                # + all_summary["probe_errors"]
            )
        # when iupac enable, after degenerated, some of the different site will become the same
        final_all_summary = all_summary.groupby(
            [
                "fwd_site",
                "fwd_errors",
                "rev_site",
                "rev_errors",
                "total_errors",
                "total_count_seqdb",
                "total_amplicon",
                "total_amplion_of_seqdb_pct",
            ]
        ).sum()

        final_all_summary.sort_values(
            ["amplicon_detected"], axis=0, ascending=[False], inplace=True
        )
        final_all_summary = final_all_summary.reset_index()
        # print(type(final_all_summary))
        # print(final_all_summary)

        # final_all_summary['amplicon_detected_of_total_amplicon_pct']
        final_all_summary["amplicon_detected_of_total_amplicon_pct"] = round(
            final_all_summary["amplicon_detected_of_total_amplicon_pct"], 2
        )
        final_all_summary["total_amplion_of_seqdb_pct"] = round(
            final_all_summary["total_amplion_of_seqdb_pct"], 2
        )
        final_all_summary.to_csv(
            path_to_assay_report,
            sep="\t",
            index=True,
            # index_label="variant_id",
        )



def main():
    """
    Main entry point for your project.
    Args:
        args : list
            A of arguments as if they were input in the command line. Leave it
            None to use sys.argv.
    """

    parser = get_parser()
    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:

        parser.print_help()
        parser.exit()


if __name__ == "__main__":
    main()
