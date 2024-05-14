import argparse
import os
from Bio import SeqIO

## Updated August 8, 2023 by Heather Ward
##    Add an optional prefix argument; used for naming complete contig files


def get_args():
    """
    Get arguments from command line with argparse.
    """
    parser = argparse.ArgumentParser(
        prog="Fasta-Make-Long-Seq-Bins.py",
        description="""Write fasta seqs to individual fasta files if longer than length threshhold.""",
    )

    parser.add_argument(
        "-p", "--prefix", required=False, help="Output file name prefix"
    )
    parser.add_argument(
        "-i", "--input_fasta", required=True, help="The input fasta file."
    )
    parser.add_argument(
        "-l",
        "--length",
        required=False,
        type=int,
        default=500000,
        help="Minimum length of contigs.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Name of output directory to write individual fastas.",
    )
    parser.add_argument(
        "-b",
        "--bins_contigs",
        required=True,
        help="Name of output file with bin number/contig names.",
    )

    return parser.parse_args()


def make_outdir(outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)


def parse_lengths_write_fastas(fasta, length, outdir, bins_contigs, prefix):
    print("Writing bins for contigs > {:,} bp...".format(length))
    count_label = 0
    contig_bin_dict = {}
    for rec in SeqIO.parse(fasta, "fasta"):
        if len(rec) >= length:
            count_label += 1
            outname = os.path.join(
                outdir, "{}complete.{}.fa".format(prefix, count_label)
            )
            contig_bin_dict["{}complete.{}".format(prefix, count_label)] = rec.id
            with open(outname, "a") as fh:
                fh.write(rec.format("fasta"))
    print("\tDone. Wrote {} contigs.".format(count_label))
    print("Writing bin:contig file...")
    with open(bins_contigs, "a") as fh:
        for k, v in contig_bin_dict.items():
            fh.write("{}\t{}\n".format(k, v))
    print("\tDone.")


def main():
    args = get_args()
    make_outdir(args.outdir)

    prefix = args.prefix + "." if args.prefix else ""

    parse_lengths_write_fastas(
        args.input_fasta, args.length, args.outdir, args.bins_contigs, prefix
    )


if __name__ == "__main__":
    main()
