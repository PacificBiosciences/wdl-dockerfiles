#!/usr/bin/env python3

from argparse import ArgumentParser
import json


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
            "id": cohort_info["cohort_id"],
            "phenotypes": cohort_info["phenotypes"],
            "affecteds": affected_samples,
            "unaffecteds": unaffected_samples,
        }
    ]

    with open(output_yaml_file, "w+") as f:
        yaml.dump(parsed_data, f)


def parse_trio(cohort_json_file):
    """
    Parse trio information from a cohort.
    Validates that i) all samples in a trio are present in the cohort; ii) exactly one trio is found in the cohort.
    Writes the index of the child, father, and mother samples to child_, father_, and mother_index.txt.
    Args:
        cohort_json_file (str): Path to the cohort JSON file
    """
    cohort_info = _load_cohort_info(cohort_json_file)

    samples_in_cohort = {
        sample["sample_id"]: sample_index
        for sample_index, sample in enumerate(cohort_info["samples"])
    }

    if len(samples_in_cohort) > 3:
        print(
            "[WARN] Found more than three samples in the cohort; attempting to extract a single valid trio"
        )

    children_with_parents = []
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
            children_with_parents.append(
                {
                    "child": {
                        "sample_id": child_id,
                        "sample_index": samples_in_cohort[child_id],
                    },
                    "father": {
                        "sample_id": father_id,
                        "sample_index": samples_in_cohort[father_id],
                    },
                    "mother": {
                        "sample_id": mother_id,
                        "sample_index": samples_in_cohort[mother_id],
                    },
                }
            )

    if len(children_with_parents) == 0:
        raise SystemExit(
            "[ERROR] Expected a trio with child and both parent samples present in the cohort; missing at least one sample in the trio"
        )
    elif len(children_with_parents) > 1:
        raise SystemExit(
            f"[ERROR] Expected a single trio, but found more than one valid trio in the cohort\n{children_with_parents}"
        )
    else:
        trio_sample = children_with_parents[0]
        print(
            f"[INFO] Found a valid trio for child sample {trio_sample['child']['sample_id']}; father sample {trio_sample['father']['sample_id']}; mother sample {trio_sample['mother']['sample_id']}"
        )
        child_index = trio_sample["child"]["sample_index"]
        father_index = trio_sample["father"]["sample_index"]
        mother_index = trio_sample["mother"]["sample_index"]

        with open("child_index.txt", "w") as f:
            f.write(f"{child_index}\n")
            print(f"[INFO] Wrote child index [{child_index}] to file child_index.txt")
        with open("father_index.txt", "w") as f:
            f.write(f"{father_index}\n")
            print(
                f"[INFO] Wrote father index [{father_index}] to file father_index.txt"
            )
        with open("mother_index.txt", "w") as f:
            f.write(f"{mother_index}\n")
            print(
                f"[INFO] Wrote mother index [{mother_index}] to file mother_index.txt"
            )


def main(args):
    if args.output_yaml_file:
        json_to_yaml(args.cohort_json, args.output_yaml_file)
    elif args.parse_trio:
        parse_trio(args.cohort_json)


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
        "--parse_trio",
        action="store_true",
        help="Validate that the cohort contains a complete trio, and output child [child_index.txt], father [father_index.txt], and mother [mother_index.txt] indices for the cohort",
    )
    options.add_argument(
        "-y",
        "--write_cohort_yaml",
        dest="output_yaml_file",
        type=str,
        help="Write the cohort information as a yaml file",
    )

    args = parser.parse_args()

    main(args)
