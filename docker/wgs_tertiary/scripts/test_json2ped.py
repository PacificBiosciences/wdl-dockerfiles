#!/usr/bin/env python3

from json2ped import parse_sample, parse_family


def test_parse_sample():
    f_id = "f"
    s = {"sample_id": "s", "sex": "MALE", "affected": True}
    assert parse_sample(f_id, s) == ["f", "s", ".", ".", "1", "2"]


def test_parse_family_trio():
    f = {
        "family_id": "f",
        "samples": [
            {"sample_id": "s", "sex": "MALE", "father_id": "d", "mother_id": "m", "affected": True},
            {"sample_id": "d", "sex": "MALE", "affected": False},
            {"sample_id": "m", "sex": "FEMALE", "affected": False}
        ],
    }
    assert parse_family(f) == [
        ["f", "s", "d", "m", "1", "2"],
        ["f", "d", ".", ".", "1", "1"],
        ["f", "m", ".", ".", "2", "1"],
    ]


def test_sex_unknown():
    f_id = "f"
    s = {"sample_id": "s", "affected": True}
    assert parse_sample(f_id, s)[4] == "."


def test_status_unknown():
    f_id = "f"
    s = {"sample_id": "s", "sex": "MALE"}
    assert parse_sample(f_id, s)[5] == "0"


def test_sex_strings():
    """Reported sex.
    (1=male; 2=female; .=unknown)"""
    f = {
        "family_id": "f",
        "samples": [
            {"sample_id": "a", "sex": "MALE",},
            {"sample_id": "a", "sex": "Male",},
            {"sample_id": "b", "sex": "M",},
            {"sample_id": "b", "sex": "m",},
            {"sample_id": "c", "sex": "FEMALE", },
            {"sample_id": "c", "sex": "Female",},
            {"sample_id": "d", "sex": "F",},
            {"sample_id": "d", "sex": "f",},
            {"sample_id": "e", "sex": "UNKNOWN",},
            {"sample_id": "f",}
        ],
    }
    assert [_[4] for _ in parse_family(f)] == ["1", "1", "1", "1", "2", "2", "2", "2", ".", "."]


def test_status_strings():
    """Field 6: phenotype (1=unaffected; 2=affected, 0=missing)"""
    f = {
        "family_id": "f",
        "samples": [
            {"sample_id": "a","affected": True,},
            {"sample_id": "b","affected": False,},
            {"sample_id": "b"},
        ],
    }
    assert [_[5] for _ in parse_family(f)] == ["2", "1", "0"]
