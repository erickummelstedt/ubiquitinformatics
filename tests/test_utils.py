import pytest
import json
import copy
from copy import deepcopy
import sys

# home_dir = os.path.expanduser('~')
local_path = '/home/erickummelstedt/lecodebase/ubiquitinformatics'
sys.path.insert(0, local_path)

from src.utils import \
    match_assertion_error_contains,\
    all_strings_exist, \
    all_strings_exist_in_list

def test_match_all_parts_present():
    msg = "Missing sites: K48. Invalid format: 1234"
    parts = ["Missing sites", "K48", "Invalid format"]
    assert match_assertion_error_contains(msg, parts)

def test_match_some_parts_missing():
    msg = "Only this is present"
    parts = ["Only", "not present"]
    assert not match_assertion_error_contains(msg, parts)

@pytest.mark.parametrize("substrings, error_string, expected", [
    # ✅ 1. All substrings present
    (["hello", "world"], "hello world", True),

    # ❌ 2. Some substrings missing
    (["hello", "python"], "hello world", False),

    # ⚠️ 3. No substrings at all (empty list)
    ([], "hello world", True),

    # ❌ 4. Target string is empty
    (["hello"], "", False),

    # ❌ 5. Case-sensitive mismatch
    (["Hello", "World"], "hello world", False),

    # ✅ 6. Full match on exact substrings
    (["quick", "brown", "fox"], "The quick brown fox jumps", True),

    # ❌ 7. One invalid substring
    (["quick", "browns"], "The quick brown fox", False),

    # ✅ 8. Overlapping substrings
    (["own", "fo"], "The quick brown fox", True),

    # ✅ 9. Repeated substrings
    (["fox", "fox"], "The fox is clever like a fox", True),

    # ✅ 10. Special characters
    (["$", "%"], "Discounts up to $50 and 10% off!", True),

    # ✅ 11. Multiline strings
    (["Line", "two", "three"], "Line one\nLine two\nLine three", True),
])
def test_all_strings_exist(substrings, error_string, expected):
    assert all_strings_exist(substrings, error_string) == expected


@pytest.mark.parametrize("error_strings, expected_substrings, expected", [
    # ✅ All matches
    (
        ["apple banana", "orange mango", "grape fruit"],
        [["apple", "banana"], ["orange", "mango"], ["grape", "fruit"]],
        [True, True, True]
    ),

    # ❌ One group fails (missing substring)
    (
        ["apple banana", "orange mango", "grape fruit"],
        [["apple", "banana"], ["orange"], ["grape", "kiwi"]],
        [True, True, False]
    ),

    # ❌ Empty group
    (
        ["apple banana", "orange mango"],
        [["apple"], []],
        [True, True]  # all() on empty iterable returns True
    ),

    # ❌ Empty target string
    (
        ["", "orange mango"],
        [["apple"], ["orange"]],
        [False, True]
    ),

    # ❌ Case sensitivity
    (
        ["Apple Banana", "ORANGE Mango"],
        [["apple"], ["orange"]],
        [False, False]
    ),

    # ✅ Special characters
    (
        ["Hello, world!", "Python > Java"],
        [["Hello", ","], ["Python", ">"]],
        [True, True]
    ),

    # ❌ Substring not found
    (
        ["one two three", "four five six"],
        [["seven"], ["five", "six", "seven"]],
        [False, False]
    ),

    # ✅ Edge case: empty input
    (
        [], [],
        []
    ),
])
def test_all_strings_exist_in_list(error_strings, expected_substrings, expected):
    assert all_strings_exist_in_list(expected_substrings, error_strings) == expected


def test_length_mismatch_raises():
    error_strings = ["string one", "string two"]
    expected_substrings = [["one"], ["two"], ["extra"]]
    with pytest.raises(ValueError, match="Length of error_string and expected_substrings must be equal"):
        all_strings_exist_in_list(expected_substrings, error_strings)




