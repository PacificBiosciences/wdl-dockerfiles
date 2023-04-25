#!/usr/bin/env python3

import csv
from operator import itemgetter
from argparse import ArgumentParser

quickviewColumns = [
    "SAMPLE",
    "PCT_OFF_BAIT",
    "MEAN_BAIT_COVERAGE",
    "MEAN_TARGET_COVERAGE",
    "FOLD_ENRICHMENT",
    "ZERO_CVG_TARGETS_PCT",
    "PCT_EXC_DUPE",
    "FOLD_80_BASE_PENALTY",
    "PCT_TARGET_BASES_20X",
    "PCT_TARGET_BASES_30X",
    "AT_DROPOUT",
    "GC_DROPOUT",
]


def main(args):
    quickgetter = itemgetter(*quickviewColumns)
    with open(args.consolidated_hs_metrics_tsv, newline="") as hsmetrics, open(
        args.output_quickview_tsv, "w", newline=""
    ) as quickview:
        reader = csv.DictReader(hsmetrics, dialect="unix", delimiter="\t")
        writer = csv.DictWriter(
            quickview, quickviewColumns, dialect="unix", delimiter="\t", quoting=0
        )
        writer.writeheader()
        for row in reader:
            writer.writerow(dict(zip(quickviewColumns, quickgetter(row))))


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Create a quickview TSV from a consolidated batch HS metrics TSV"
    )

    parser.add_argument(
        "-b",
        "--consolidated-hs-metrics-tsv",
        type=str,
        help="Consolidated batch HS metrics TSV",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-quickview-tsv",
        type=str,
        help="Ouptut file to write consolidated quickview TSV to",
        required=True,
    )

    args = parser.parse_args()
    main(args)
