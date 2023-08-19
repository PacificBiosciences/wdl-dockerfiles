#!/usr/bin/env python3
"""
Split uBAM file into chunks of chunksize reads
"""

__version__ = "0.1.0"


import argparse
from os.path import basename

import pysam


def main(args):
    save = pysam.set_verbosity(0)  # suppress [E::idx_find_and_load]
    with pysam.AlignmentFile(args.ubam, check_sq=False) as bamin:
        pysam.set_verbosity(save)  # restore warnings
        recordcount = 0
        chunk = 0
        for record in bamin:
            if recordcount == 0:
                bamout = pysam.AlignmentFile(
                    basename(args.ubam).replace(".bam", f".chunk_{chunk}.bam"),
                    "wb",
                    template=bamin,
                )
            bamout.write(record)
            recordcount += 1
            if recordcount == int(args.chunksize):
                bamout.close()
                recordcount = 0
                chunk += 1
    bamout.close()


if __name__ == "__main__":
    """This is executed when run from the command line"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ubam", help="Unaligned reads BAM")
    parser.add_argument("chunksize", help="Maximum number of reads per chunk")
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    args = parser.parse_args()
    main(args)
