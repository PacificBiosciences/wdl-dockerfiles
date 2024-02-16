#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from argparse import ArgumentParser


def _parse_args():
    parser = ArgumentParser(description="Plot read categories for a batch")

    parser.add_argument(
        "-r",
        "--batch-read-data-csv",
        type=str,
        help="Read data CSV containing read data for all samples in the batch",
        required=True,
    )
    parser.add_argument(
        "-l",
        "--lima-demultiplex-report",
        type=str,
        help="Demultiplex report output by lima",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--batch-duplicate-lengths-csv",
        type=str,
        help="Duplicate read lengths csv for the batch",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        type=str,
        help="Output directory to save files into",
        default=".",
    )

    args = parser.parse_args()
    return args


def main():
    try:
        readCsv = snakemake.input.readCsv
        lima = snakemake.input.lima
        dups = snakemake.input.dups
        outdir = snakemake.params.odir
    except NameError:
        args = _parse_args()
        readCsv = args.batch_read_data_csv
        lima = args.lima_demultiplex_report
        dups = args.batch_duplicate_lengths_csv
        outdir = args.out_dir

    allDemuxReads = pd.read_csv(readCsv)
    limaReport = pd.read_csv(lima, sep="\t", usecols=["ReadLengths", "PassedFilters"])
    dups = pd.read_csv(dups, usecols=["sample", "length", "category"])

    # load data
    onTarget = allDemuxReads.query('target != "."').assign(category="OnTarget Unique")[
        ["sample", "category", "length"]
    ]

    offTarget = allDemuxReads.query('target == "." and ismapped==1').assign(
        category="OffTarget Unique"
    )[["sample", "category", "length"]]

    unmapped = allDemuxReads.query("ismapped==0").assign(category="Unmapped Unique")[
        ["sample", "category", "length"]
    ]

    noDemux = (
        limaReport.query("PassedFilters==0")
        .assign(category="Failed Demux")
        .drop("PassedFilters", axis=1)
        .rename(columns={"ReadLengths": "length"})
    )

    # merge
    allData = pd.concat(
        [onTarget, offTarget, unmapped, dups, noDemux], ignore_index=True
    )

    # plot
    hue_order = [
        "OnTarget Unique",
        "OffTarget Unique",
        "Unmapped Unique",
        "Duplicates",
        "Failed Demux",
    ]

    f, (ax0, ax1) = plt.subplots(
        1, 2, figsize=(10, 8), gridspec_kw={"width_ratios": [1, 5]}
    )

    # histogram
    sns.histplot(
        data=allData,
        x="length",
        hue_order=hue_order,
        hue="category",
        multiple="stack",
        binrange=(0, 15000),
        ax=ax1,
    )

    # side plot
    data = (
        allData.groupby("category")
        .size()
        .reindex(hue_order[::-1])
        .rename("reads")
        .reset_index()
    )
    data["yieldFrac"] = data.reads.transform(lambda x: x / x.sum()).fillna(0)

    bottom = 0
    colors = sns.color_palette()[:5][::-1]

    for i, row in data.iterrows():
        ax0.bar(1, row.yieldFrac, bottom=bottom, color=colors[i], alpha=0.75)
        ax0.annotate(
            f"{row.yieldFrac:.2f}",
            (1, bottom + 0.5 * row.yieldFrac),
            ha="center",
            va="center",
        )
        bottom += row.yieldFrac

    ax0.set_xticks([])
    ax0.set_ylabel("Total Fraction")

    plt.tight_layout()

    f.savefig(f"{outdir}/read_categories.png", dpi=400)

    # write readlength legnth stats table
    allData.fillna("noSampleData").groupby(
        ["category", "sample"]
    ).length.describe().dropna().astype(int).to_csv(
        f"{outdir}/read_length_by_sample.csv"
    )


if __name__ == "__main__":
    main()
