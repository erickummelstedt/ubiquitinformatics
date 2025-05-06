import pytest
import json
import copy
from copy import deepcopy
import sys
import logging

# home_dir = os.path.expanduser('~')
# local_path = '/home/erickummelstedt/lecodebase/ubiquitinformatics/src/main.py'
local_path = '/Users/ekummelstedt/le_code_base/ubiquitinformatics/back_end'
sys.path.insert(0, local_path)

from src.utils import \
    match_assertion_error_contains,\
    all_strings_exist, \
    all_strings_exist_in_list, \
    convert_json_to_dict, \
    inject_fasta_sequence_at_chain,\
    inject_protein_key,\
    inject_branching_sites

from tests.test_data import \
    five_level_nested_ubiquitin_

from src.logging_utils import \
    log_branching_details,\
    log_end_of_branching, \
    log_protein_details, \
    log_end_of_protein

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




@pytest.mark.parametrize("input_data, expected_output", [
    # ✅ Case 1: Input already dictionary
    (
        {"protein": "1ubq", "chain_number": 1},
        {"protein": "1ubq", "chain_number": 1}
    ),
    # ✅ Case 2: Valid JSON string with double quotes
    (
        '{"protein": "1ubq", "chain_number": 1}',
        {"protein": "1ubq", "chain_number": 1}
    ),
    # ✅ Case 3: Valid JSON string with single quotes (handled)
    (
        "{'protein': '1ubq', 'chain_number': 1}",
        {"protein": "1ubq", "chain_number": 1}
    ),
    # ✅ Case 4: JSON string with spaces and formatting
    (
        '{ "protein" : "1ubq" , "chain_number" : 1 }',
        {"protein": "1ubq", "chain_number": 1}
    )
])
def test_convert_json_to_dict_valid(input_data, expected_output):
    result = convert_json_to_dict(input_data)
    assert result == expected_output, f"Expected {expected_output}, got {result}"


@pytest.mark.parametrize("invalid_input", [
    "invalid_json",                  # Not JSON
    "{'protein': '1ubq', 'chain_number': }",  # Malformed JSON
    123,                            # Integer input
    ["not", "a", "dict"],            # List input
    None                             # None input
])
def test_convert_json_to_dict_invalid(invalid_input):
    if isinstance(invalid_input, str):
        with pytest.raises(ValueError, match="Invalid JSON format"):
            convert_json_to_dict(invalid_input)
    else:
        with pytest.raises(TypeError, match="Input must be a dictionary or a JSON string"):
            convert_json_to_dict(invalid_input)

@pytest.mark.parametrize("branch, working_dictionary, expected_logs", [

    # ✅ Test 1: Standard branch with normal lysine
    (
        {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
        {"chain_number": 1},
        [
            "===== START OF LYSINE SITE =====",
            "Chain Number: 1",
            "Lysine Site: K48",
            "Sequence ID: FAG(K)QLE",
            "Lysine Conjugation: "
        ]
    ),

    # ✅ Test 2: Branch with protecting group SMAC
    (
        {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
        {"chain_number": 2},
        [
            "===== START OF LYSINE SITE =====",
            "Chain Number: 2",
            "Lysine Site: K29",
            "Sequence ID: VKA(K)IQD",
            "Lysine Conjugation: SMAC"
        ]
    ),

    # ✅ Test 3: Branch with protecting group ABOC
    (
        {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": "ABOC"},
        {"chain_number": 3},
        [
            "===== START OF LYSINE SITE =====",
            "Chain Number: 3",
            "Lysine Site: K6",
            "Sequence ID: IFV(K)TLT",
            "Lysine Conjugation: ABOC"
        ]
    ),

    # ✅ Test 4: Branch with special characters
    (
        {"site_name": "K🔥", "sequence_id": "🔥SEQ🔥", "children": "ABOC"},
        {"chain_number": 4},
        [
            "===== START OF LYSINE SITE =====",
            "Chain Number: 4",
            "Lysine Site: K🔥",
            "Sequence ID: 🔥SEQ🔥",
            "Lysine Conjugation: ABOC"
        ]
    ),

    # ✅ Test 5: Large chain number stress test
    (
        {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
        {"chain_number": 99999},
        [
            "===== START OF LYSINE SITE =====",
            "Chain Number: 99999",
            "Lysine Site: K48",
            "Sequence ID: FAG(K)QLE",
            "Lysine Conjugation: "
        ]
    ),

])
def test_log_branching_details(caplog, branch, working_dictionary, expected_logs):
    """Test logging output of log_branching_details with various branches."""
    with caplog.at_level(logging.INFO):
        log_branching_details(branch, working_dictionary, context=None)  # context is not used

    logs = [record.message for record in caplog.records]

    for expected_line in expected_logs:
        assert any(expected_line in log for log in logs), f"Expected log line '{expected_line}' not found in logs: {logs}"


# ✅ Test 6: Branch missing 'children' field
def test_log_branching_details_missing_children(caplog):
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE"}  # Missing 'children'
    working_dictionary = {"chain_number": 1}

    with pytest.raises(KeyError, match='children'):
        log_branching_details(branch, working_dictionary, context=None)
    

# ✅ Test 7: Branch missing 'sequence_id' field
def test_log_branching_details_missing_sequence_id(caplog):
    branch = {"site_name": "K48", "children": "SMAC"}  # Missing 'sequence_id'
    working_dictionary = {"chain_number": 1}

    with pytest.raises(KeyError, match='sequence_id'):
        log_branching_details(branch, working_dictionary, context=None)
    

# ✅ Test 8: Branch missing 'site_name' field
def test_log_branching_details_missing_site_name(caplog):
    branch = {"sequence_id": "FAG(K)QLE", "children": "ABOC"}  # Missing 'site_name'
    working_dictionary = {"chain_number": 1}

    with pytest.raises(KeyError, match='site_name'):
        log_branching_details(branch, working_dictionary, context=None)

def test_log_end_of_branching(caplog):
    """
    Test that `log_end_of_branching` logs the correct message.
    """
    with caplog.at_level(logging.INFO):
        log_end_of_branching()

    # Get all log messages
    logs = [record.message for record in caplog.records]

    # Check that the expected message is in the logs
    assert "===== END OF LYSINE SITE =====" in logs, f"Expected log message not found. Got: {logs}"

@pytest.fixture
def sample_working_dictionary():
    return {
        "protein": "1ubq",
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""}
        ]
    }

@pytest.fixture
def sample_context():
    return {
        "chain_length_list": [76],
        "chain_number_list": [1]
    }

def test_log_protein_details(caplog, sample_working_dictionary, sample_context):
    """
    Test that `log_protein_details` logs correct protein information.
    """
    with caplog.at_level(logging.INFO):
        log_protein_details(sample_working_dictionary, sample_context)

    logs = [record.message for record in caplog.records]

    # Check for expected messages
    expected_logs = [
        " ===== START OF PROTEIN ===== ",
        f"Protein: {sample_working_dictionary['protein']}",
        f"Sequence: {sample_working_dictionary['FASTA_sequence']}",
        f"Chain Length List: {sample_context['chain_length_list']}",
        f"Chain Length: {sample_working_dictionary['chain_length']}",
        f"Chain Number List: {sample_context['chain_number_list']}",
        f"Chain Number: {sample_working_dictionary['chain_number']}",
        f"Branching Sites: {sample_working_dictionary.get('branching_sites', [])}"
    ]

    for expected in expected_logs:
        assert expected in logs, f"Expected log message '{expected}' not found in logs."

def test_log_end_of_protein(caplog, sample_working_dictionary):
    """
    Test that `log_end_of_protein` logs correct end message.
    """
    with caplog.at_level(logging.INFO):
        log_end_of_protein(sample_working_dictionary)

    logs = [record.message for record in caplog.records]
    expected_log = f"===== END OF PROTEIN - UBI NUMBER: {sample_working_dictionary['chain_number']} ====="

    assert expected_log in logs, f"Expected log message '{expected_log}' not found in logs."


@pytest.mark.parametrize("chain_number, new_fasta_sequence", [
    (2, "MODSEQ_CHAIN_2"),
    (3, "MODSEQ_CHAIN_3"),
    (4, "MODSEQ_CHAIN_4"),
    (5, "MODSEQ_CHAIN_5")
])
def test_inject_fasta_sequence_at_chain(chain_number, new_fasta_sequence):
    """
    Test that inject_fasta_sequence_at_chain correctly updates the FASTA_sequence
    at the target chain number.
    """
    nested_ub = copy.deepcopy(five_level_nested_ubiquitin_)

    inject_fasta_sequence_at_chain(
        nested_ub["branching_sites"],
        target_chain_number=chain_number,
        new_fasta_sequence=new_fasta_sequence
    )

    def get_fasta_at_chain(branches, target, current=1):
        """
        Recursively find the FASTA_sequence at a specified chain number.
        """
        for site in branches:
            if isinstance(site.get("children"), dict):
                if site["children"].get("chain_number") == target:
                    return site["children"]["FASTA_sequence"]
                result = get_fasta_at_chain(site["children"]["branching_sites"], target, current + 1)
                if result:
                    return result
        return None

    modified_fasta = get_fasta_at_chain(nested_ub["branching_sites"], chain_number)

    assert modified_fasta == new_fasta_sequence, (
        f"Expected FASTA sequence at chain {chain_number} to be '{new_fasta_sequence}', got '{modified_fasta}'"
    )


@pytest.mark.parametrize("chain_number, key, value", [
    (2, "new_key", "new_value"),
    (3, "chain_label", 777),
    (5, "note", "terminal protein"),
])
def test_inject_protein_key_add_key(chain_number, key, value):
    """Ensure a key is correctly injected at a specific chain number."""
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    inject_protein_key(nested["branching_sites"], target_chain_number=chain_number, key=key, value=value)

    def get_key_value(branches, target):
        for site in branches:
            child = site.get("children")
            if isinstance(child, dict):
                if child.get("chain_number") == target:
                    return child.get(key)
                result = get_key_value(child.get("branching_sites", []), target)
                if result:
                    return result
        return None

    assert get_key_value(nested["branching_sites"], chain_number) == value


@pytest.mark.parametrize("chain_number, key", [
    (2, "chain_length"),
    (5, "FASTA_sequence"),
    (4, "protein")
])
def test_inject_protein_key_remove_key(chain_number, key):
    """Ensure a key is correctly removed from a specific chain (non-recursive per case)."""
    nested = copy.deepcopy(five_level_nested_ubiquitin_)

    # Apply removal
    inject_protein_key(nested["branching_sites"], target_chain_number=chain_number, key=key, remove=True)

    # Manual path lookup for each chain level
    if chain_number == 2:
        # Chain 2 is at K48 of root (chain 1)
        target = nested["branching_sites"][6]["children"]
    elif chain_number == 4:
        # Chain 4 is at K6 of chain 3 → K63 of chain 2 → K48 of root
        target = (
            nested["branching_sites"][6]["children"]  # chain 2
            ["branching_sites"][7]["children"]        # chain 3
            ["branching_sites"][0]["children"]        # chain 4
        )
    elif chain_number == 5:
        # Chain 5 is at K11 of chain 3 → K63 of chain 2 → K48 of root
        target = (
            nested["branching_sites"][6]["children"]  # chain 2
            ["branching_sites"][7]["children"]        # chain 3
            ["branching_sites"][2]["children"]        # chain 5
        )
    else:
        raise ValueError(f"Unsupported chain number for test: {chain_number}")

    assert key not in target, f"Key '{key}' was not removed from chain {chain_number}"


def test_inject_protein_key_stops_after_first_match():
    """Ensure it modifies only the first matching chain (no duplicates in tree)."""
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    inject_protein_key(nested["branching_sites"], 3, "unique_marker", "FOUND")

    def count_key_occurrences(branches, key_name):
        count = 0
        for site in branches:
            child = site.get("children")
            if isinstance(child, dict):
                if key_name in child:
                    count += 1
                count += count_key_occurrences(child.get("branching_sites", []), key_name)
        return count

    assert count_key_occurrences(nested["branching_sites"], "unique_marker") == 1


def test_inject_protein_key_no_crash_on_missing_children():
    """Ensure no crash if children are not present or not a dict."""
    structure = {
        "branching_sites": [
            {"site_name": "K48", "children": ""},  # no child dict
            {"site_name": "K63", "children": "SMAC"}  # protected group, not dict
        ]
    }

    # Should not raise even though there are no valid children
    inject_protein_key(structure["branching_sites"], 2, "key", "value")

# ===================================
# TESTS: inject_branching_sites
# This is a simplified version of the five_level_nested_ubiquitin_ structure
# ===================================

# Branching sites we'll inject for testing
custom_sites = [
    {"site_name": "K99", "sequence_id": "ZZZ(K)ZZZ", "children": ""},
    {"site_name": "K100", "sequence_id": "YYY(K)YYY", "children": ""}
]

def test_inject_chain_2_branching_sites():
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    inject_branching_sites(nested["branching_sites"], 2, custom_sites)

    # Navigate manually to chain 2
    chain2 = nested["branching_sites"][6]["children"]
    assert chain2["chain_number"] == 2
    assert chain2["branching_sites"] == custom_sites

def test_inject_chain_3_branching_sites():
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    chain2 = nested["branching_sites"][6]["children"]
    chain3 = chain2["branching_sites"][7]["children"]
    inject_branching_sites(nested["branching_sites"], 3, custom_sites)

    assert chain3["chain_number"] == 3
    assert chain3["branching_sites"] == custom_sites

def test_inject_chain_4_branching_sites():
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    chain2 = nested["branching_sites"][6]["children"]
    chain3 = chain2["branching_sites"][7]["children"]
    chain4 = chain3["branching_sites"][1]["children"]
    inject_branching_sites(nested["branching_sites"], 4, custom_sites)

    assert chain4["chain_number"] == 4
    assert chain4["branching_sites"] == custom_sites

def test_inject_chain_5_branching_sites():
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    chain2 = nested["branching_sites"][6]["children"]
    chain3 = chain2["branching_sites"][7]["children"]
    chain5 = chain3["branching_sites"][2]["children"]
    inject_branching_sites(nested["branching_sites"], 5, custom_sites)

    assert chain5["chain_number"] == 5
    assert chain5["branching_sites"] == custom_sites

def test_no_injection_when_chain_does_not_exist():
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    original = copy.deepcopy(nested)

    inject_branching_sites(nested["branching_sites"], 99, custom_sites)

    # Expect structure unchanged since chain 99 doesn't exist
    assert nested == original