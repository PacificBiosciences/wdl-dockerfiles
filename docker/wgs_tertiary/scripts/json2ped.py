#!/usr/bin/env python3
"""
Convert Family JSON structure to tab-delimited PLINK pedigree (PED) format.

Output PED columns:
1. family_id
2. sample_id
3. father_id (. for unknown)
4. mother_id (. for unknown)
5. sex (1=male; 2=female; .=unknown)
6. phenotype (1=unaffected; 2=affected)
"""

__version__ = "0.1.0"

import json
import csv
import sys


SEX = {"MALE": "1", "M": "1", "FEMALE": "2", "F": "2"}
STATUS = {False: "1", True: "2"}


def parse_sample(family_id, sample):
    """For a sample struct, return a list of PED fields."""
    return [
        family_id,
        sample["sample_id"],
        sample.get("father_id", "."),
        sample.get("mother_id", "."),
        SEX.get(sample.get("sex", ".").upper(), "."),  # all cases accepted
        STATUS.get(sample.get("affected"), "0"),
    ]


def parse_family(family):
    """For a family struct, return a list of lists of PED fields for each sample."""
    family_id = family["family_id"]
    samples = []
    for sample in family["samples"]:
        samples.append(parse_sample(family_id, sample))
    return samples


def write_ped(samples):
    """Write PED format to stdout."""
    tsv_writer = csv.writer(sys.stdout, delimiter="\t")
    for sample in samples:
        tsv_writer.writerow(sample)


def main():
    with open(sys.argv[1], "r") as family:
        samples = parse_family(json.load(family))
        write_ped(samples)


if __name__ == "__main__":
    if sys.argv[1] in ["-v", "--version"]:
        print(__version__)
        sys.exit(0)
    main()
