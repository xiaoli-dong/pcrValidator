import os
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError

__version__ = "0.1.0"
__author__ = "Xiaoli Dong"


def get_parser():

    # Disable default help
    parser = argparse.ArgumentParser(
        add_help=False,
        epilog="For example:\n python3 ./bin/fetch_gbk.py -q txid64320[Organism:noexp]  -d nucleotide -e test@example.com -r gb",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        "--query",
        "-q",
        type=str,
        required=True,
        help="user defined query string to retrieve data from NCBI database",
    )
    required.add_argument(
        "--dbname",
        "-d",
        type=str,
        required=True,
        default="nucleotide",
        help="NCBI database name",
    )
    required.add_argument(
        "-email",
        "-e",
        type=str,
        required=True,
        help="user email to access NCBI database",
    )
    # prog give the program basename without path but __file__ gives the full path
    version = parser.prog + " " + __version__
    optional.add_argument("--version", "-v", action="version", version=version)

    optional.add_argument(
        "--prefix",
        "-p",
        type=str,
        default="download",
        help="prefix of the output file names",
    )
    optional.add_argument(
        "--rettype",
        "-r",
        type=str,
        default="fasta",
        choices=["fasta", "gb"],
        help="retrieval type when for ncbi esearch ",
    )
    optional.add_argument(
        "--retmax",
        "-m",
        type=int,
        default="500",
        help="Total number of records from the input set to be retrieved from ncbi every time min > 0 and max < 100,000",
    )

    return parser


def download(query, dbname, email, rettype, retmax, prefix):

    Entrez.email = email
    search_handle = Entrez.esearch(db=dbname, term=query, usehistory="y")
    search_results = Entrez.read(search_handle)
    print(search_results)
    search_handle.close()
    count = int(search_results["Count"])
    print("Found %i results" % count)
    if os.path.exists(prefix + "." + rettype):
        print(prefix + "." + rettype + " file already exist, will not download ...")
        return

    out_handle = open(prefix + "." + rettype, "w")
    for start in range(0, count, retmax):
        end = min(count, start + retmax)
        print("Going to download record %i to %i" % (start + 1, end))
        fetch_handle = Entrez.efetch(
            db=dbname,
            rettype=rettype,
            retmode="text",
            retstart=start,
            retmax=retmax,
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"],
        )

        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)
    out_handle.close()


def extract_fasta_and_metadata(input, rettype, fasta_seq_out, source_feature_out):
    count = 0
    with open(fasta_seq_out, "w") as fasta_out, open(
        source_feature_out, "w"
    ) as meta_out:

        for seq_record in SeqIO.parse(input, rettype):
            id = seq_record.id
            if seq_record.seq:
                try:
                    # print(seq_record.id)
                    # print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
                    # print(seq_record.seq)
                    r = SeqIO.write(seq_record, fasta_out, "fasta")
                    for feat in seq_record.features:
                        if feat.type == "source":
                            meta_out.write(f"{id}")
                            for k, v in feat.qualifiers.items():
                                # v is a list
                                # print(f"{id} /{k}={' '.join(v)}", end="")
                                meta_out.write(f" /{k}={' '.join(v)}")
                                # print("xxxxxxxxxxxxxxxxx")
                            # print("\n")
                            meta_out.write("\n")
                    count += r

                except UndefinedSequenceError:
                    print(seq_record.id + " sequence content is undefined")
                except:
                    print("Error while writing sequence:  " + seq_record.id)

    print(f"In total write {count} sequences\n")


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
    download(
        args.query, args.dbname, args.email, args.rettype, args.retmax, args.prefix
    )
    if args.rettype != "fasta":
        extract_fasta_and_metadata(
            args.prefix + "." + args.rettype,
            args.rettype,
            args.prefix + ".fasta",
            args.prefix + "_metadata.csv",
        )


if __name__ == "__main__":
    main()
