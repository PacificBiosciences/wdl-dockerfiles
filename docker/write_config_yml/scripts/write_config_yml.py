#!/usr/bin/env python3

from argparse import ArgumentParser
import json
import yaml


def parse_json(json_file):
    """
    Parse json data to the required key/value format
    Args:
        json_file (str): Path to the cohort.json file
    Returns:
        (dict) cohort information dictionary
    """
    with open(json_file, "r") as f:
        data = json.load(f)

    affected_samples = list()
    unaffected_samples = list()
    for sample in data["samples"]:
        parents = [
            parent
            for parent in [sample["father_id"], sample["mother_id"]]
            if parent is not None
        ]
        sex = sample["sex"]
        if sex not in ["MALE", "FEMALE"]:
            raise SystemExit(f"Invalid sex [{sex}]; must be one of ['MALE', 'FEMALE']")
        sample_info = {
            "id": sample["sample_id"],
            "sex": sample["sex"],
            "parents": parents,
        }
        if sample["affected"] == True:
            affected_samples.append(sample_info)
        elif sample["affected"] == False:
            unaffected_samples.append(sample_info)

    parsed_data = [
        {
            "id": data["cohort_id"],
            "phenotypes": data["phenotypes"],
            "affecteds": affected_samples,
            "unaffecteds": unaffected_samples,
        }
    ]
    return parsed_data


def main(args):
    json_dict = parse_json(args.cohort_json)

    with open(args.output_yml, "w+") as f:
        yaml.dump(json_dict, f)


if __name__ == "__main__":
    parser = ArgumentParser(description="Write a cohort configuration yml file")

    parser.add_argument(
        "-j", "--cohort_json", type=str, help="Cohort JSON file", required=True
    )
    parser.add_argument(
        "-y", "--output_yml", type=str, help="Output yml file to write", required=True
    )

    args = parser.parse_args()
    main(args)
