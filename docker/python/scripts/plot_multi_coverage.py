#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from math import ceil
from argparse import ArgumentParser


def _parse_args():
    parser = ArgumentParser(description="Plot batch coverage")

    parser.add_argument(
        "-t",
        "--targets-coverage-summary-csv",
        type=str,
        help="Coverage CSV summary for all samples in batch across target regions bed",
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
        targets_coverage_summary_csv = snakemake.input[0]
        out_dir = snakemake.params.odir
    except NameError:
        args = _parse_args()
        targets_coverage_summary_csv = args.targets_coverage_summary_csv
        out_dir = args.out_dir

    DPI = 400

    data = pd.read_csv(targets_coverage_summary_csv)
    # replace "." with "off-target"
    ot = "off-target"
    data.target = data.target.str.replace(".", ot, regex=False)

    ####
    # multi-sample coverage per target
    ####

    order = sorted(data.target.unique())
    ncols = ceil(len(order) ** 0.5)

    g = sns.FacetGrid(
        data=data,
        sharex=False,
        col="target",
        col_wrap=ncols,
        col_order=order,
        hue="sample",
    )

    g.map(plt.plot, "start", "coverage").set_xlabels("chr start pos").set_ylabels(
        "Coverage"
    ).add_legend().savefig(f"{out_dir}/multi_coverage_by_target.png", dpi=DPI)


if __name__ == "__main__":
    main()
