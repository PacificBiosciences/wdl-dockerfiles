#!/usr/bin/env python3
"""
Measure the consistency of non-reference kmers between
all pairs of provided kmer counts.  Typically each kmer
count file is from a SMRT Cell.  The kmers presents in
different SMRT Cells from the same sample should be consistent
with each other.  It is preferred that kmer counts be filtered
to only report modimers.
"""

__version__ = "0.2.0"


import argparse
import gzip
import io

THRESHOLD = 0.03  # empirically determined threshold for kmer consistency


def read_kmers(kmers_tsv, solid_count=5, return_solid=False):
    """Read kmers and counts from a tsv file.
    Return a set of all kmers and optionally a set of solid kmers.
    """
    kmers = set()
    solid_kmers = set()
    with io.TextIOWrapper(gzip.open(kmers_tsv)) as f:
        for line in f:
            kmer, count = line.rstrip("\n").split()
            count = int(count)
            kmers.add(kmer)
            # a kmer is solid if it is seen at least solid_count times
            if count >= solid_count:
                solid_kmers.add(kmer)
    if return_solid:
        return (kmers, solid_kmers)
    else:
        return (kmers, )


def kmer_inconsistency(ds1_kmers, ds2_kmers, refkmers):
    """
    For each pair of datasets, count the number of shared and unique
    reference and non-reference kmers.  Use the count of unique reference
    kmers to adjust for undersampling (i.e. low coverage) since they are
    likely shared between any two humans.
    """
    # unpack kmers tuples
    ds1_allkmers, ds1_solidkmers = ds1_kmers
    ds2_allkmers, ds2_solidkmers = ds2_kmers

    # kmers solid in at least one dataset
    solidkmers = ds1_solidkmers | ds2_solidkmers
    # "shared" kmers solid in at least one dataset and present in both
    sharedsolidkmers = solidkmers & ds1_allkmers & ds2_allkmers
    # "unique" kmers solid in one dataset and absent in the other
    uniquekmers = solidkmers - sharedsolidkmers
    # all "shared" kmers whether or not they are solid
    sharedkmers = ds1_allkmers & ds2_allkmers

    # count separately for reference and non-reference kmers
    ref_shared = len(sharedkmers & refkmers)
    ref_unique = len(uniquekmers & refkmers)
    refkmer_inconsistency = (
        0 if ref_unique == 0 else ref_unique / (ref_shared + ref_unique)
    )
    nonref_shared = len(sharedkmers - refkmers)
    nonref_unique = len(uniquekmers - refkmers)
    nonrefkmer_inconsistency = (
        0 if nonref_unique == 0 else nonref_unique / (nonref_shared + nonref_unique)
    )

    # output adjusted kmer consistency (nonref - ref)
    return max(0, nonrefkmer_inconsistency - refkmer_inconsistency)


def main(args):
    # Read the reference kmers/modimers
    refkmers = read_kmers(args.ref_kmers_tsv, return_solid=False)[0]

    # Read the kmers for each sample
    datasetkmers = list()
    for ds in args.dataset_kmers_tsv:
        datasetkmers.append(read_kmers(ds, return_solid=True))

    print("movieA\tmovieB\tadjusted_nonref_inconsistency\tconsistent")
    for ds1ix in range(len(args.dataset_kmers_tsv)):
        movie1 = (".").join(args.dataset_kmers_tsv[ds1ix].split("/")[-1].split(".")[0:-3])
        for ds2ix in range(ds1ix + 1, len(args.dataset_kmers_tsv)):
            movie2 = (".").join(args.dataset_kmers_tsv[ds2ix].split("/")[-1].split(".")[0:-3])

            adjusted_nonrefkmer_inconsistency = kmer_inconsistency(
                datasetkmers[ds1ix], datasetkmers[ds2ix], refkmers
            )
            inconsistent = (
                "YES" if adjusted_nonrefkmer_inconsistency < THRESHOLD else "NO"
            )

            print(
                f"{movie1}\t{movie2}\t{adjusted_nonrefkmer_inconsistency:0.5f}\t{inconsistent}"
            )


if __name__ == "__main__":
    """This is executed when run from the command line"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "ref_kmers_tsv", help="Reference genome kmer counts (kmer<TAB>count)"
    )
    parser.add_argument(
        "dataset_kmers_tsv", nargs="+", help="Kmer counts (kmer<TAB>count)"
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    args = parser.parse_args()
    main(args)