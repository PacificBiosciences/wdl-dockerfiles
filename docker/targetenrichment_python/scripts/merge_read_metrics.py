#!/usr/bin/env python3

import pandas as pd
from argparse import ArgumentParser


def _parse_args():
    parser = ArgumentParser(description="Merge read metrics")

    parser.add_argument("-s", "--sample-id", type=str, help="Sample ID", required=True)
    parser.add_argument(
        "-r",
        "--read-metrics-csv",
        type=str,
        help="Path to sample reads CSV",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--read-metrics-targets-csv",
        type=str,
        help="Path to read targets CSV",
        required=True,
    )
    parser.add_argument(
        "-e",
        "--exons-per-read-csv",
        type=str,
        help="Path to exons per read CSV",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-merged-read-metrics",
        type=str,
        help="Path to output merged read metrics CSV",
        required=True,
    )

    args = parser.parse_args()
    return args


def main():
    try:
        csvs = snakemake.input
        sample_id = snakemake.params.sample
        output_merged_read_metrics = snakemake.output[0]
    except NameError:
        args = _parse_args()
        sample_id = args.sample_id
        csvs = [
            args.read_metrics_csv,
            args.read_metrics_targets_csv,
            args.exons_per_read_csv,
        ]
        output_merged_read_metrics = args.output_merged_read_metrics

    res = pd.concat(
        [pd.read_csv(csv, index_col="readname") for csv in csvs], axis=1
    ).reset_index()
    res.insert(0, "sample", sample_id)
    res.to_csv(output_merged_read_metrics, index=False)


if __name__ == "__main__":
    main()
