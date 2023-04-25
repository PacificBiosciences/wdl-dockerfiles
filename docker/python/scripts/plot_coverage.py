#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def _parse_args():
    parser = ArgumentParser(description="Plot coverage data for a sample")

    parser.add_argument(
        "-t",
        "--targets-base-coverage-csv",
        type=str,
        help="Coverage CSV for sample across target regions bed",
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
        targets_base_coverage_csv = snakemake.input[0]
        out_dir = snakemake.params.odir
    except NameError:
        args = _parse_args()
        targets_base_coverage_csv = args.targets_base_coverage_csv
        out_dir = args.out_dir

    DPI = 200

    coverage = pd.read_csv(targets_base_coverage_csv)

    order = sorted(coverage.target.unique())
    maxCols = int(len(order) ** 0.5) + 1

    plt.figure(figsize=(40, 40))

    g = sns.FacetGrid(
        data=coverage,
        sharex=False,
        col="target",
        col_wrap=maxCols,
        col_order=order,
        hue="target",
    )

    g.map(plt.plot, "start", "coverage")

    g.set_xlabels("chr start pos")
    g.set_ylabels("Coverage")

    g.savefig(f"{out_dir}/coverage_by_target.png", dpi=DPI)


if __name__ == "__main__":
    main()
