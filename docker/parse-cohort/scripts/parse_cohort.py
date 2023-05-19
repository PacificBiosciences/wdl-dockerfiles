#!/usr/bin/env python3

from argparse import ArgumentParser
import json

__version__ = "1.0.2"


def _load_cohort_info(cohort_json_file):
    """
    Load a cohort json file
    Args:
        cohort_json_file (str): Path to the cohort JSON file
    Returns:
        (dict): A dict of cohort information
    """
    with open(cohort_json_file, "r") as f:
        cohort_info = json.load(f)

    return cohort_info


def json_to_yaml(cohort_json_file, output_yaml_file):
    """
    Parse a cohort info dict to the required key/value format, and write the output to a yaml file
    Args:
        cohort_json_file (str): Path to the cohort JSON file
        output_yaml_file (str): File to write yaml information to
    """
    import yaml

    cohort_info = _load_cohort_info(cohort_json_file)

    affected_samples = list()
    unaffected_samples = list()
    for sample in cohort_info["samples"]:
        parents = [
            parent
            for parent in [sample["father_id"], sample["mother_id"]]
            if parent is not None
        ]
        sex = sample.get("sex")
        if sex not in ["MALE", "FEMALE", None]:
            raise SystemExit(
                f"Invalid sex [{sex}]; must be one of ['MALE', 'FEMALE', null]"
            )
        sample_info = {
            "id": sample["sample_id"],
            "parents": parents,
        }
        if sex:
            sample_info["sex"] = sex

        if sample["affected"]:
            affected_samples.append(sample_info)
        else:
            unaffected_samples.append(sample_info)

    parsed_data = [
        {
            "id": cohort_info["cohort_id"],
            "phenotypes": cohort_info["phenotypes"],
            "affecteds": affected_samples,
            "unaffecteds": unaffected_samples,
        }
    ]

    with open(output_yaml_file, "w+") as f:
        yaml.dump(parsed_data, f)


def parse_families(cohort_json_file):
    """
    Parse family/trio information from a cohort.
    Validates that all samples in a trio are present in the cohort
    Writes the children, father, and mother indices for each valid trio/family to families.json
    Args:
        cohort_json_file (str): Path to the cohort JSON file
    """
    cohort_info = _load_cohort_info(cohort_json_file)

    samples_in_cohort = {
        sample["sample_id"]: sample_index
        for sample_index, sample in enumerate(cohort_info["samples"])
    }

    trios = dict()
    for sample in cohort_info["samples"]:
        child_id = sample["sample_id"]
        father_id = sample["father_id"]
        mother_id = sample["mother_id"]

        # Check that both parental IDs are defined, and both parents are found in this cohort
        if (
            father_id
            and mother_id
            and set(samples_in_cohort.keys()).issuperset(set([father_id, mother_id]))
        ):
            child_index = samples_in_cohort[child_id]

            parent_key = f"{father_id}_{mother_id}"
            if parent_key in trios:
                trios[parent_key]["child_indices"].append(child_index)
            else:
                father_index = samples_in_cohort[father_id]
                mother_index = samples_in_cohort[mother_id]
                trios[parent_key] = {
                    "child_indices": [child_index],
                    "father_index": father_index,
                    "mother_index": mother_index,
                }

    if len(trios) == 0:
        raise SystemExit(
            "[ERROR] Expected at least one trio with child and both parent samples present in the cohort; missing at least one sample in the trio"
        )

    else:
        with open("families.json", "w+") as f:
            json.dump(list(trios.values()), f)


def main(args):
    if args.output_yaml_file:
        json_to_yaml(args.cohort_json, args.output_yaml_file)
    elif args.parse_families:
        parse_families(args.cohort_json)


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Parse a cohort to validate and retrieve trio sample indices or to write cohort information to a yaml file"
    )

    parser.add_argument(
        "-j", "--cohort_json", type=str, help="Cohort JSON file", required=True
    )

    options = parser.add_mutually_exclusive_group(required=True)
    options.add_argument(
        "-t",
        "--parse_families",
        action="store_true",
        help="Validate that the cohort contains at minimum one complete trio, and output sample indices for each valid trio in the cohort as families.json (array of valid families)",
    )
    options.add_argument(
        "-y",
        "--write_cohort_yaml",
        dest="output_yaml_file",
        type=str,
        help="Write the cohort information as a yaml file",
    )
    options.add_argument(
        "--version",
        action="version",
        version=f"parse_cohort.py version {__version__}"
    )

    args = parser.parse_args()

    main(args)
