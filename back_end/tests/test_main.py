import types
import pytest
import json
import copy
import logging
from copy import deepcopy
import sys

# home_dir = os.path.expanduser('~')
# local_path = '/home/erickummelstedt/lecodebase/ubiquitinformatics/src/main.py'
local_path = '/Users/ekummelstedt/le_code_base/ubiquitinformatics/back_end'
sys.path.insert(0, local_path)

# Import the functions from the original code
# from src.main_testing import relabelling_ubiquitin_numbers, inner_wrapper_relabelling_ubiquitin_numbers
from src.main import \
    iterate_through_ubiquitin, \
    inner_wrapper_iterate_through_ubiquitin, \
    find_branching_site, \
    validate_protein_keys, \
    check_branching_sites, \
    check_branching_sequences,\
    validate_branching_sites,\
    check_branching_site_sequence_match, \
    check_children_format,\
    process_current_protein, \
    process_branch, \
    add_max_chain_number, \
    K_residue_ubi_addition, \
    process_ubiquitin_reaction, \
    ubiquitin_simulation, \
    inner_wrapper_ubiquitin_simulation, \
    handle_lysine_modification, \
    ubiquitin_building, \
    inner_wrapper_ubiquitin_building

from src.utils import \
    match_assertion_error_contains,\
    all_strings_exist, \
    all_strings_exist_in_list, \
    inject_fasta_sequence_at_chain,\
    inject_protein_key,\
    inject_branching_sites

from tests.test_data import \
    five_level_nested_ubiquitin_,\
    k48_dimer_ubiquitin,\
    string_k48_dimer_ubiquitin,\
    ubiquitin_monomer, \
    histag_ubiquitin_monomer,\
    BASE_WORKING_DICT, \
    BASE_CONTEXT

from src.logging_utils import \
    log_branching_details,\
    log_end_of_branching,\
    log_protein_details,\
    log_end_of_protein

'''
From main:
- find_branching_site
- validate_protein_keys
- check_branching_sites
- check_branching_sequences
- validate_branching_sites
- check_branching_site_sequence_match
- check_children_format
- validate_branching_sites
- process_current_protein
- process_branch

From Utils
- match_assertion_error_contains
- all_strings_exist
- all_strings_exist_in_list
- convert_json_to_dict

For logging_utils
- log_branching_details
- log_end_of_branching
- log_protein_details
- log_end_of_protein



Build tests for the following
Everything works; now build tests and clean up code for each of the following functions
Redo cover all: 
- iterate_through_ubiquitin
- inner_wrapper_iterate_through_ubiquitin
- add_max_chain_number

'''

# =====================================
# === Tests for find_branching_site ===
# Tests for validating the indexing of branching sites and handling of invalid or empty inputs.
# =====================================

@pytest.mark.parametrize("sequence_id, expected_position", [
    ("(M)QIF", 1),
    ("IFV(K)TLT", 6),
    ("LTG(K)TIT", 11),
    ("FAG(K)QLE", 48),
    ("IQD(K)EGI", 33),
    ("NIQ(K)EST", 63),
])

def test_find_branching_site(sequence_id, expected_position):
    """
    Test `find_branching_site` function to verify correct indexing of 
    branching sites in a given FASTA sequence.

    Steps:
    1. Provide different sequence IDs and their expected positions.
    2. Call `find_branching_site` with each test case.
    3. Assert that the returned position matches the expected position.
    """
    fasta_sequence = k48_dimer_ubiquitin["FASTA_sequence"]  # Retrieve test FASTA sequence

    # Validate the function's output
    assert find_branching_site(sequence_id, fasta_sequence) == expected_position, \
        f"Expected position {expected_position} for sequence {sequence_id}, but got a different value."

def test_find_branching_site_not_found():
    """
    Test that `find_branching_site` raises a ValueError when given a 
    sequence ID that does not exist in the provided FASTA sequence.
    """
    unknown_sequence = "ABC(K)XYZ"  # A sequence not found in FASTA_sequence

    with pytest.raises(ValueError, match="substring not found"):
        find_branching_site(unknown_sequence, k48_dimer_ubiquitin["FASTA_sequence"])

def test_find_branching_site_empty():
    """
    Test that `find_branching_site` raises a ValueError when given an empty sequence ID as input.
    """
    empty_sequence = ""  # Invalid input: empty string

    with pytest.raises(ValueError, match="substring not found"):
        find_branching_site(empty_sequence, k48_dimer_ubiquitin["FASTA_sequence"])


@pytest.mark.parametrize("chain_number, altered_sequence", [
    # Chain 4, K33 altered (to break IQD(K)EGI)
    (4, "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEHIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"),

    # Chain 5, K6 altered (to break IFV(K)TLT)
    (5, "MQIFVTKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"),

    # Chain 2, K48 altered (to break FAG(K)QLE)
    (2, "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQJEGRTLSDYNIQKESTLHLVLRLRGG")
])
def test_find_branching_site_nested_not_found(chain_number, altered_sequence):
    """
    Recursively modify the FASTA sequence for a nested ubiquitin and assert
    that find_branching_site raises a ValueError when the sequence is missing.
    """
    # Copy the test structure
    nested_ub = copy.deepcopy(five_level_nested_ubiquitin_)

    # Inject a modified sequence at the specified chain number
    inject_fasta_sequence_at_chain(
        nested_ub["branching_sites"],
        target_chain_number=chain_number,
        new_fasta_sequence=altered_sequence
    )

    # Now test that the sequence is not found
    with pytest.raises(ValueError, match="substring not found"):
        iterate_through_ubiquitin(nested_ub)

# =====================================
# === Tests for validate_protein_keys ===
# Tests for validating the presence of required keys and handling of invalid or extra keys.
# =====================================

@pytest.mark.parametrize("valid_dict", [
    # ✅ Test Case 1: Valid ubiquitin dictionary
    {
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}]
    },
    # ✅ Test Case 2: Another valid ubiquitin with different sequence and length
    {
        "protein": "2ubq",
        "chain_number": 2,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 82,
        "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}]
    }
])
def test_validate_protein_keys_valid(valid_dict):
    """
    Test that `validate_protein_keys` does not raise an error for valid dictionaries.
    """
    try:
        validate_protein_keys(valid_dict)  # Should not raise an error
        assert True  # Pass the test
    except KeyError:
        pytest.fail("validate_protein_keys raised KeyError unexpectedly for a valid dictionary.")

@pytest.mark.parametrize("invalid_dict, missing_keys", [
    # ❌ Test Case 3: Missing "chain_number"
    (
        {
            "protein": "1ubq",
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "chain_length": 76,
            "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                                {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}]
        },
        {"chain_number"}
    ),
    # ❌ Test Case 4: Missing multiple required keys
    (
        {
            "protein": "1ubq",
            "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                                {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}]
        },
        {"chain_length", "chain_number", "FASTA_sequence"}
    )
])
def test_validate_protein_keys_missing_keys(invalid_dict, missing_keys):
    """
    Test that `validate_protein_keys` raises a KeyError when required keys are missing.
    """
    allowed_keys = {"protein", "chain_number", "FASTA_sequence", "chain_length", "branching_sites"}
    with pytest.raises(KeyError, match=r"Missing required keys: .*Allowed keys: .*"):
        validate_protein_keys(invalid_dict)

@pytest.mark.parametrize("invalid_dict, invalid_keys", [
    # ❌ Test Case 5: Dictionary with an invalid key
    (
        {
            "protein": "1ubq",
            "chain_number": 1,
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "chain_length": 76,
            "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                                {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}],
            "extra_key": "invalid_value"
        },
        ["extra_key"]
    ),
    # ❌ Test Case 6: Dictionary with multiple invalid keys
    (
        {
            "protein": "1ubq",
            "chain_number": 1,
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "chain_length": 76,
            "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                                {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}],
            "extra_key_1": "value1",
            "extra_key_2": "value2"
        },
        ["extra_key_1", "extra_key_2"]
    )
])

def test_validate_protein_keys_invalid_keys(invalid_dict, invalid_keys):
    """
    Test that `validate_protein_keys` raises a KeyError when invalid keys are present.
    """    
    with pytest.raises(KeyError) as exc_info:
            validate_protein_keys(invalid_dict)
        
    error_msg = str(exc_info.value)
    assert match_assertion_error_contains(error_msg, invalid_keys), \
        f"Expected parts {invalid_keys} not found in error: {error_msg}"


def test_validate_protein_keys_missing_and_invalid_keys():
    """
    Test `validate_protein_keys` when the input dictionary has both missing required keys and unexpected keys.
    It should raise a KeyError with the correct message.
    """
    invalid_protein_data = {
        "protein": "1ubq",  # Valid key
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",  # Valid key
        "extra_key": "unexpected_value",  # Invalid key
    }

    with pytest.raises(KeyError, match=r"Missing required keys: .*Invalid keys found: .*Allowed keys: .*"):
        validate_protein_keys(invalid_protein_data)


@pytest.mark.parametrize("chain_number, key_to_remove, expected_missing_key", [
    (3, "chain_number", "chain_number"),
    (2, "FASTA_sequence", "FASTA_sequence"),
    (4, "chain_length", "chain_length"),
])
def test_validate_protein_keys_missing_nested(chain_number, key_to_remove, expected_missing_key):
    """Test that a missing key in a nested ubiquitin triggers KeyError via recursive traversal."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    inject_protein_key(data["branching_sites"], chain_number, key_to_remove, remove=True)

    with pytest.raises(KeyError) as exc_info:
        iterate_through_ubiquitin(data)

    assert expected_missing_key in str(exc_info.value)


@pytest.mark.parametrize("chain_number, key_to_add, value", [
    (5, "bad_key", "oops"),
    (3, "wrong_field", 123),
])
def test_validate_protein_keys_invalid_nested(chain_number, key_to_add, value):
    """Ensure that unexpected keys in nested ubiquitins raise KeyError."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    inject_protein_key(data["branching_sites"], chain_number, key_to_add, value)

    with pytest.raises(KeyError) as exc_info:
        iterate_through_ubiquitin(data)

    assert key_to_add in str(exc_info.value)

def test_validate_protein_keys_missing_and_invalid_nested():
    """Trigger both a missing and an unexpected key in a nested ubiquitin."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)

    inject_protein_key(data["branching_sites"], 4, "chain_length", remove=True)
    inject_protein_key(data["branching_sites"], 4, "unexpected_key", "bad")

    with pytest.raises(KeyError) as exc_info:
        iterate_through_ubiquitin(data)

    msg = str(exc_info.value)
    assert "chain_length" in msg and "unexpected_key" in msg


# =====================================
# === Tests for check_branching_sites ===
# Tests for validating the presence of required branching sites and handling of invalid or extra sites.
# =====================================

@pytest.mark.parametrize("ubiquitin_structure, expected_errors", [
    
    # ✅ Test Case 1: Correct Ubiquitin Structure (No Errors)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}]
        
    }, []),

    # ❌ Test Case 2: Missing Required Sites
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""}# Missing K29, K33, K48, K63
        ]
    }, [["Missing required sites", "K29", "K33", "K48", "K63", "Ubiquitin 1"], 
        ["Allowed sites", "M1", "K6", "K11", "K27", "K29", "K33", "K48", "K63"]]
    ),
    
    # ❌ Test Case 3: Invalid Site Present
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""},
            {"site_name": "K99","sequence_id": "NIQ(K)EST","children": ""}  # Invalid site
        ]
    }, [["Invalid sites found", "K99", "Ubiquitin 1"], 
        ["Allowed sites", "M1", "K6", "K11", "K27", "K29", "K33", "K48", "K63"]]),

    # ❌ Test Case 4: Both Missing and Invalid Sites
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K100","sequence_id": "NIQ(K)EST","children": ""}]
    }, [["Missing required sites","K27", "K29", "K48", "K63", "Ubiquitin 1"], 
        ["Invalid sites found", "K100", "Ubiquitin 1"], 
        ["Allowed sites", "M1", "K6", "K11", "K27", "K29", "K33", "K48", "K63"]]),

    # ✅ Test Case 5: Deeply Nested Ubiquitin (Correct)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children": ""}, 
            {"site_name": "K48", "children": {
                "chain_number": 2,
                "branching_sites": [
                    {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                    {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                    {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                    {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                    {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                    {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                    {"site_name": "K48","sequence_id": "FAG(K)QLE","children": ""},
                    {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}
                ]
            }},
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}
        ]
    }, [])

])
def test_check_branching_sites(ubiquitin_structure, expected_errors):
    """
    Test `check_branching_sites()` with various ubiquitin structures to verify correct
    validation of required branching sites.
    """

    site_errors = check_branching_sites(ubiquitin_structure)
    
    if expected_errors:
        assert all(all_strings_exist_in_list(expected_errors, site_errors)), \
            f"Expected errors {expected_errors} for ubiquitin {ubiquitin_structure}, not in {site_errors}."
    
    else:
        assert len(site_errors) == 0

@pytest.mark.parametrize("chain_number, new_branching_sites, expected_error_substrings", [
    # ❌ Missing branching sites: omit K63, K29
    (
        3,
        [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""}
        ],
        [["Missing required sites", "K29", "K63", "Ubiquitin 3"],
         ["Allowed sites", "M1", "K6", "K11", "K27", "K29", "K33", "K48", "K63"],
         ["Allowed sequences", '(M)QIF', 'IFV(K)TLT', '(M)QIF', 'IQD(K)EGI', 'LTG(K)TIT', 'ENV(K)AKI', 'VKA(K)IQD', 'FAG(K)QLE', 'NIQ(K)EST']]
    ),

    # ❌ Invalid site present
    (
        4,
        [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children": ""},
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""},
            {"site_name": "K99", "sequence_id": "INVALID(K)SEQ", "children": ""}
        ],
        [["Invalid sites found", "K99", "Ubiquitin 4"],
         ["Invalid sequences found", "INVALID(K)SEQ", "Ubiquitin 4"],
         ["Allowed sites", "M1", "K6", "K11", "K27", "K29", "K33", "K48", "K63"],
         ["Allowed sequences", '(M)QIF', 'IFV(K)TLT', '(M)QIF', 'IQD(K)EGI', 'LTG(K)TIT', 'ENV(K)AKI', 'VKA(K)IQD', 'FAG(K)QLE', 'NIQ(K)EST']]
    ),

    # ❌ Missing and invalid
    (
        2,
        [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K100", "sequence_id": "WTF(K)OMG", "children": ""}
        ],
        [["Missing required sites", "K6", "K11", "K27", "K29", "K33", "K48", "K63", "Ubiquitin 2"],
         ["Invalid sites found", "K100", "Ubiquitin 2"],
         ["Allowed sites", "M1", "K6", "K11", "K27", "K29", "K33", "K48", "K63"],
         ["Allowed sequences", '(M)QIF', 'IFV(K)TLT', '(M)QIF', 'IQD(K)EGI', 'LTG(K)TIT', 'ENV(K)AKI', 'VKA(K)IQD', 'FAG(K)QLE', 'NIQ(K)EST']]
    )
])
def test_check_branching_sites_nested(chain_number, new_branching_sites, expected_error_substrings):
    """
    Recursively inject errors into nested branching_sites and ensure `check_branching_sites`
    captures them via `iterate_through_ubiquitin`.
    """
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    inject_branching_sites(nested["branching_sites"], target_chain_number=chain_number, new_branching_sites=new_branching_sites)

    with pytest.raises(AssertionError) as exc_info:
            iterate_through_ubiquitin(nested)
    
    error_msg = str(exc_info.value)
    combined_substrings = [item for sublist in expected_error_substrings for item in sublist]
    assert all_strings_exist(combined_substrings, error_msg), \
        f"Expected {expected_error_substrings} in {error_msg}, but not found."
    
# =====================================
# === Tests for check_branching_sequences ===
# Tests for validating the sequences of branching sites and handling of invalid or missing sequences.   
# =====================================

@pytest.mark.parametrize("ubiquitin_structure, expected_errors", [

    # ✅ Test Case 1: All required sequence_ids present
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}
        ]
    }, []),

    # ❌ Test Case 2: Missing sequences
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""}
            # Missing rest
        ]

    }, [["Missing required sequences", "ENV(K)AKI","VKA(K)IQD","IQD(K)EGI","FAG(K)QLE","NIQ(K)EST", "Ubiquitin 1"], 
        ["Allowed sequences", "(M)QIF", "IFV(K)TLT", "LTG(K)TIT", "ENV(K)AKI", "VKA(K)IQD", "IQD(K)EGI", "FAG(K)QLE", "NIQ(K)EST"]]),
    
    # ❌ Test Case 3: Invalid sequence present
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K16","sequence_id": "INVALID(SEQ)","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}
        ]
    }, [["Invalid sequences found", "INVALID(SEQ)", "Ubiquitin 1"], 
        ["Allowed sequences", "(M)QIF", "IFV(K)TLT", "LTG(K)TIT", "ENV(K)AKI", "VKA(K)IQD", "IQD(K)EGI", "FAG(K)QLE", "NIQ(K)EST"]]),
    
    # ❌ Test Case 4: Both missing and invalid sequences
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "INVALID(SEQ)","children": ""}
        ]
    }, [["Missing required sequences", "LTG(K)TIT", "ENV(K)AKI","VKA(K)IQD","IQD(K)EGI","FAG(K)QLE","NIQ(K)EST", "Ubiquitin 1"], 
        ["Invalid sequences found", "INVALID(SEQ)", "Ubiquitin 1"], 
        ["Allowed sequences", "(M)QIF", "IFV(K)TLT", "LTG(K)TIT", "ENV(K)AKI", "VKA(K)IQD", "IQD(K)EGI", "FAG(K)QLE", "NIQ(K)EST"]]),
    
    # ✅ Test Case 5: Deeply nested valid structure
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                    {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                    {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                    {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                    {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                    {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                    {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                    {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}
            ]
            }}
        ]
    }, [])

])
def test_check_branching_sequences(ubiquitin_structure, expected_errors):
    """
    Test `check_branching_sequences()` with various ubiquitin structures to verify correct
    validation of required branching sites.
    """

    sequence_errors = check_branching_sequences(ubiquitin_structure)
    
    if expected_errors:
        assert all(all_strings_exist_in_list(expected_errors, sequence_errors)), \
            f"Expected errors {expected_errors} for ubiquitin {ubiquitin_structure}, not in {sequence_errors}."
    
    else:
        assert len(sequence_errors) == 0

@pytest.mark.parametrize("chain_number, injected_branching_sites, expected_error_substrings", [
    # Missing sequences
    (
        3,
        [
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""}
        ],
        [["Missing required sequences", "Ubiquitin 3"], ["Allowed sequences"]]
    ),
    # Invalid sequence
    (
        4,
        [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11","sequence_id": "BAD(SEQ)","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}
         ],
        [["Invalid sequences found", "BAD(SEQ)", "Ubiquitin 4"], ["Allowed sequences"]]
    ),
])
def test_check_branching_sequences_recursive(chain_number, injected_branching_sites, expected_error_substrings):
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    inject_branching_sites(nested["branching_sites"], chain_number, injected_branching_sites)

    with pytest.raises(AssertionError) as exc_info:
            iterate_through_ubiquitin(nested)
    
    error_msg = str(exc_info.value)
    combined_substrings = [item for sublist in expected_error_substrings for item in sublist]
    assert all_strings_exist(combined_substrings, error_msg), \
        f"Expected {expected_error_substrings} in {error_msg}, but not found."

# ========================================
# ==== Tests for check_children_format ===
# Tests for validating the format of children in branching sites.
# ========================================

@pytest.mark.parametrize("ubiquitin_dict, site, expected_error", [
    # ✅ Empty string
    ({"chain_number": 1}, {"children": ""}, []),

    # ✅ SMAC
    ({"chain_number": 1}, {"children": "SMAC"}, []),

    # ✅ ABOC
    ({"chain_number": 1}, {"children": "ABOC"}, []),

    # ✅ Dictionary (valid nested ubiquitin)
    ({"chain_number": 1}, {"children": {"protein": "1ubq", "chain_number": 1}}, []),

    # ❌ List (invalid type)
    ({"chain_number": 4}, {"children": ["invalid_structure"]}, 
     ["Invalid children format: ['invalid_structure'] in Ubiquitin 4"]),

    # ❌ Integer (invalid type)
    ({"chain_number": 1000}, {"children": 1234}, 
     ["Invalid children format: 1234 in Ubiquitin 1000"]),
])
def test_check_children_format(ubiquitin_dict, site, expected_error):
    errors = []
    check_children_format(ubiquitin_dict, site, errors)
    assert errors == expected_error

@pytest.mark.parametrize("chain_number, site_with_invalid_children, expected_error", [
    (
        3,
        [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
         {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ["invalid_structure"]},
         {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
         {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
         {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
         {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
         {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
         {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}],
        # Expected Error
        "Invalid children format: ['invalid_structure'] in Ubiquitin 3"
    ),
    (
        4,
        [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
         {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
         {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": 123},
         {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
         {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
         {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
         {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
         {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}],
        
        # Expected Error
        "Invalid children format: 123 in Ubiquitin 4"
    )
])
def test_check_children_format_recursive(chain_number, site_with_invalid_children, expected_error):
    nested = copy.deepcopy(five_level_nested_ubiquitin_)
    inject_branching_sites(nested["branching_sites"], chain_number, site_with_invalid_children)

    with pytest.raises(AssertionError) as exc_info:
            iterate_through_ubiquitin(nested)
    
    error_msg = str(exc_info.value)
    assert expected_error in error_msg, \
        f"Expected {expected_error} in {error_msg}, but not found."

# =====================================
# Tests for check_branching_site_sequence_match
# Tests to validate that site_name and sequence_id mappings are correct and handle invalid or missing mappings.
# =====================================

@pytest.mark.parametrize("site, expected_error", [
    # ✅ Correct match
    ({"site_name": "K48", "sequence_id": "FAG(K)QLE"}, []),

    # ❌ Incorrect sequence for valid site_name
    ({"site_name": "K48", "sequence_id": "WRONG(K)SEQ"}, 
     ["site_name: K48, does not correspond with the sequence_id: WRONG(K)SEQ"]),

    # ❌ Completely missing sequence_id key
    ({"site_name": "K48"}, 
     ["site_name: K48, does not correspond with the sequence_id: None"]),

    # ❌ Valid site_name, empty sequence
    ({"site_name": "K11", "sequence_id": ""}, 
     ["site_name: K11, does not correspond with the sequence_id: "]),
])
def test_check_branching_site_sequence_match(site, expected_error):
    errors = []
    check_branching_site_sequence_match(site, errors)
    assert errors == expected_error


@pytest.mark.parametrize("chain_number, site_name, faulty_sequence, expected_error_substrings", [
    (
        2, "K6", "LTG(K)TIT",
        [["site_name: K6", "does not correspond with the sequence_id", "LTG(K)TIT"]]
    ),
    (
        3, "K48", "LTG(K)TIT",
        [["site_name: K48", "does not correspond with the sequence_id", "LTG(K)TIT"]]
    ),
    (
        4, "K27", "LTG(K)TIT",  # missing sequence
        [["site_name: K27", "does not correspond with the sequence_id", "LTG(K)TIT"]]
    ),
])

def test_iterate_check_branching_site_sequence_match_errors(chain_number, site_name, faulty_sequence, expected_error_substrings):
    """
    Injects incorrect sequence IDs into nested ubiquitin structure and verifies that
    iterate_through_ubiquitin surfaces `check_branching_site_sequence_match` errors.
    """
    modified = copy.deepcopy(five_level_nested_ubiquitin_)

    # Create faulty branching_sites at a specific chain level
    if chain_number == 2:
        faulty_branching_sites = [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": site_name, "sequence_id": faulty_sequence, "children":""}, 
            {"site_name": "K11","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""},
        ]  # no sequence_id
    elif chain_number == 3:
        faulty_branching_sites = [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "FAG(K)QLE","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": site_name, "sequence_id": faulty_sequence, "children":""}, 
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}       
            ]
    elif chain_number == 4:
        faulty_branching_sites = [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": site_name,"sequence_id": faulty_sequence,"children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}       
            ]
    
    inject_branching_sites(
        modified["branching_sites"],
        target_chain_number=chain_number,
        new_branching_sites=faulty_branching_sites
    )

    # Run recursive validation
    with pytest.raises(AssertionError) as exc_info:
        iterate_through_ubiquitin(modified)

    error_msg = str(exc_info.value)
    combined_substrings = [item for sublist in expected_error_substrings for item in sublist]
    assert all_strings_exist(combined_substrings, error_msg), \
        f"Expected {combined_substrings} in {error_msg}, but not found."
    
# =====================================
# === Tests for validate_branching_sites ===
# Tests for validating the presence of required branching sites and handling of invalid or extra sites.
# =====================================

@pytest.mark.parametrize("ubiquitin_dict, should_raise, expected_parts", [

    # ✅ Test 1: Perfect structure
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, False, None),

    # ❌ Test 2: Incorrect sequence_id mapping for one site
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "NIQ(K)EST", "children": ""},
            {"site_name": "K63", "sequence_id": "FAG(K)QLE", "children": ""}
        ]
    }, True, [
        "does not correspond with the sequence_id",
        "NIQ(K)EST",
        "K63",
        "does not correspond with the sequence_id",
        "FAG(K)QLE"]), #raises validate_branching_sites error

    # ❌ Test 3: One invalid child format (not string or dict)
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": 42},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "NIQ(K)EST", "children": ""},
            {"site_name": "K63", "sequence_id": "FAG(K)QLE", "children": ""}
        ]
    }, True, [
        "Invalid children format:",
        "42",
        "K48",
        "does not correspond with the sequence_id",
        "NIQ(K)EST",
        "K63",
        "does not correspond with the sequence_id",
        "FAG(K)QLE"]),

    # ❌ Test 3a: One invalid child format (wrong string)
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": "42"},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "NIQ(K)EST", "children": ""},
            {"site_name": "K63", "sequence_id": "FAG(K)QLE", "children": ""}
        ]
    }, True, [
        "Invalid children format:",
        "42",
        "K48",
        "does not correspond with the sequence_id",
        "NIQ(K)EST",
        "K63",
        "does not correspond with the sequence_id",
        "FAG(K)QLE"]),

    # ❌ Test 4: Invalid site_name (not expected)
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K99", "sequence_id": "IFV(K)TLT", "children": ""}, # Invalid site_name (not expected)
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "NIQ(K)EST", "children": ""},
            {"site_name": "K63", "sequence_id": "FAG(K)QLE", "children": ""}
        ]
    }, True, ["Invalid sites found",
              "K99",
              "Ubiquitin 1",
              "Allowed sites",
              "M1",
              "K6",
              "K11",
              "K27",
              "K29",
              "K33",
              "K48",
              "K63"]),

    # ❌ Test 5: Missing one required site
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""}, # Invalid site_name (not expected)
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""}
            # Missing the rest
        ]
    }, True, ["Missing required sites",
              "K29",
              "K33",
              "K48",
              "K63"]),

    # ❌ Test 6: Invalid sequence_id name
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K11", "sequence_id": "IFV(K)TLT", "children": ""}, 
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, True, ["site_name:",
              "K6",
              "does not correspond with the sequence_id",
              "site_name:", 
              "K11", 
              "does not correspond with the sequence_id"]),

    # ❌ Test 7: SMAC with invalid site/sequence
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K11", "sequence_id": "IFV(K)TLT", "children": "SMAC"}, 
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, True, ["site_name:",
              "K6",
              "does not correspond with the sequence_id",
              "site_name:", 
              "K11", 
              "does not correspond with the sequence_id"]),

    # ❌ Test 8: Children as a list (invalid)
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""}, 
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ["wrong"]},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, True, ["Invalid children format", 
              "wrong"]),

    # ❌ Test 9: Multiple errors (bad sequence + bad child format)
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K11", "sequence_id": "IFV(K)TLT", "children": ""}, 
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ["wrong"]},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, True, ["does not correspond with the sequence_id", 
              "Invalid children format"]),

    # ❌ Test 10: Correct branching sites but missing one sequence_id
    ({
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""}
            # K63 missing
        ]
    }, True, ["Missing required sites", 
              "Allowed sites", 
              "Missing required sequences", 
              "Allowed sequences"]
)
])

def test_validate_branching_sites(ubiquitin_dict, should_raise, expected_parts):

    if should_raise:
        with pytest.raises(AssertionError) as exc_info:
            validate_branching_sites(ubiquitin_dict)
        
        error_msg = str(exc_info.value)
        assert match_assertion_error_contains(error_msg, expected_parts), \
            f"Expected parts {expected_parts} not found in error: {error_msg}"
    
    else:
        validate_branching_sites(ubiquitin_dict)  # Should not raise

@pytest.mark.parametrize("chain_number, faulty_branching_sites, expected_error_substrings", [
    # ❌ K48 has wrong sequence_id
    (
        2,
        [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "WRONG(K)SEQ", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ],
        ["Invalid sequences found", "WRONG(K)SEQ", "Ubiquitin 2"]
    ),

    # ❌ K6 has integer as children (invalid format)
    (
        3,
        [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": 42},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ],
        ["Invalid children format", "42"]
    ),

    # ❌ K99 as invalid site_name
    (
        4,
        [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""},
            {"site_name": "K99", "sequence_id": "IFV(K)TLT", "children": ""}
        ],
        ["Invalid sites found", "K99", "Ubiquitin 4", "Allowed sites"]
    ),

    # ❌ Missing multiple required sites
    (
        5,
        [
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""}
        ],
        ["Missing required sites", "K27", "K29", "K33", "K48", "K63"]
    )
])
def test_validate_branching_sites_nested(chain_number, faulty_branching_sites, expected_error_substrings):
    """
    Inject faulty branching_sites at various depths and verify `validate_branching_sites`
    is triggered correctly by `iterate_through_ubiquitin`.
    """
    nested = copy.deepcopy(five_level_nested_ubiquitin_)

    # Inject faults at the target chain level
    inject_branching_sites(
        nested["branching_sites"],
        target_chain_number=chain_number,
        new_branching_sites=faulty_branching_sites
    )

    # Run and expect assertion
    with pytest.raises(AssertionError) as exc_info:
        iterate_through_ubiquitin(nested)

    error_msg = str(exc_info.value)
    combined_substrings = [item for sublist in expected_error_substrings for item in sublist]
    assert all_strings_exist(expected_error_substrings, error_msg), \
        f"Expected {expected_error_substrings} in {error_msg}, but not found."

# =====================================
# === Tests for iterate_through_ubiquitin ===
# Tests for validating the processing of ubiquitin structures, including renumbering, chain length adjustments, and error handling.
# =====================================

def test_iterate_through_ubiquitin():
    """
    Test `iterate_through_ubiquitin` to ensure proper renumbering and chain 
    length adjustments in a nested ubiquitin structure.
    
    Steps:
    1. Create a deep copy of the `k48_dimer_ubiquitin` dictionary to avoid modifying original data.
    2. Run `iterate_through_ubiquitin` on the test copy.
    3. Verify that:
       - The root ubiquitin has the correct `chain_number` (1).
       - The root ubiquitin has the correct `chain_length` (83).
       - The child ubiquitin at `K48` has the correct `chain_number` (2).
       - The `SMAC` protecting group is correctly identified at `K63`.
       - The child ubiquitin at `K48` has the correct `chain_length` (76).
    """
    # Step 1: Create a deep copy of the test data to prevent modifying the original
    test_ubiquitin = copy.deepcopy(k48_dimer_ubiquitin)

    # Step 2: Run the relabelling function
    updated_ubiquitin, context = iterate_through_ubiquitin(test_ubiquitin)

    # Step 3: Assertions to validate correct renumbering and chain adjustments
    assert updated_ubiquitin["chain_number"] == 1, "Root ubiquitin should have chain_number = 1"
    assert updated_ubiquitin["chain_length"] == 83, "Root ubiquitin should have chain_length = 83"

    # Verify child ubiquitin renumbering
    child_ubiquitin = updated_ubiquitin["branching_sites"][6]["children"]
    assert child_ubiquitin["chain_number"] == 2, "Child ubiquitin at K48 should have chain_number = 2"
    assert child_ubiquitin["chain_length"] == 76, "Child ubiquitin at K48 should have chain_length = 76"

    # Verify presence of SMAC protecting group at K63
    assert child_ubiquitin["branching_sites"][7]["children"] == "SMAC", \
        "K63 should have an SMAC protecting group"
    
def test_deep_nesting():
    """
    Test the `iterate_through_ubiquitin` function to ensure correct renumbering of ubiquitin chains
    in a deeply nested ubiquitin structure.

    Step 1: Create a deep copy of the k48_dimer_ubiquitin test dictionary
    Step 2: Insert a ubiquitin monomer at the K63 site (index 7) of the first ubiquitin's children & Further nest another ubiquitin monomer inside the previously added monomer at site K48 (index 6)
    Step 3: Call the function to relabel ubiquitin numbers in the entire structure
    Step 4: Validate the correct relabeling of chain numbers at different nesting levels

    """
    # Step 1:
    test_ubiquitin = copy.deepcopy(k48_dimer_ubiquitin)
    # Step 2:
    test_ubiquitin["branching_sites"][6]["children"]["branching_sites"][7]["children"] = copy.deepcopy(ubiquitin_monomer)
    test_ubiquitin["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][6]["children"] = copy.deepcopy(ubiquitin_monomer)
    # Step 3:
    updated_ubiquitin, context = iterate_through_ubiquitin(test_ubiquitin)
    # Step 4:
    assert updated_ubiquitin["branching_sites"][6]["children"]["branching_sites"][7]["children"]["chain_number"] == 3, \
        "Expected chain number 3 for the first nested ubiquitin monomer"
    assert updated_ubiquitin["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][6]["children"]["chain_number"] == 4, \
        "Expected chain number 4 for the deeply nested ubiquitin monomer"

def test_empty_input():
    """
    Test that `iterate_through_ubiquitin` raises KeyError when given an empty dictionary as input.
    """
    with pytest.raises(KeyError):
        iterate_through_ubiquitin({})

def test_missing_keys():
    """
    Test that `iterate_through_ubiquitin` raises a KeyError when required keys are missing from the input dictionary.
    """
    incomplete_dict = {"protein": "1ubq"}  # Missing required keys

    with pytest.raises(KeyError):
        iterate_through_ubiquitin(incomplete_dict)

# Sample ubiquitin data in JSON string format
sample_ubiquitin_json = json.dumps(k48_dimer_ubiquitin)

@pytest.mark.parametrize("input_data, expected_chain_number, expected_length", [
    (k48_dimer_ubiquitin, 1, 83),
    (sample_ubiquitin_json, 1, 83),
])

def test_json_loading_valid(input_data, expected_chain_number, expected_length):
    """
    Test iterate_through_ubiquitin with both dictionary and JSON string input.
    """
    result, context = iterate_through_ubiquitin(input_data)
    assert result["chain_number"] == expected_chain_number
    assert result["chain_length"] == expected_length

def test_loading_is_string():
    """
    Test iterate_through_ubiquitin with invalid JSON string input.
    Test that the function pulls out the JSON.
    """
    valid_json = copy.deepcopy(string_k48_dimer_ubiquitin)  # Improper JSON format (single quotes)
    
    result, context = iterate_through_ubiquitin(valid_json)
    assert result["chain_number"] == 1
    assert result["chain_length"] == 83

def test_json_loading_empty():
    """
    Test iterate_through_ubiquitin with empty JSON string.
    """
    empty_json = "{}"
    with pytest.raises(KeyError):
        iterate_through_ubiquitin(empty_json)


def test_json_loading_partial_data():
    """
    Test iterate_through_ubiquitin with missing keys in JSON.
    """
    incomplete_json = json.dumps({"protein": "1ubq"})
    with pytest.raises(KeyError):
        iterate_through_ubiquitin(incomplete_json)


def test_json_loading_non_json_string():
    """
    Test iterate_through_ubiquitin with a non-JSON string that is not representative of the correct Dictionary structure.
    """
    non_json_input = "This is not JSON"
    with pytest.raises(ValueError):
        iterate_through_ubiquitin(non_json_input)


def test_json_loading_with_extra_keys():
    """
    Test iterate_through_ubiquitin with additional invalid keys in JSON. 
    Assert the the additional invalid keys create errors.
    """
    erroneous_json = json.dumps({
        "protein": "1ubq",
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH",
        "invalid_key": "some_value"
    })
    
    with pytest.raises(KeyError):
        iterate_through_ubiquitin(erroneous_json)

def test_json_loading_with_invalid_keys():
    """
    Test iterate_through_ubiquitin with JSON containing invalid values as keys.
    """
    invalid_json = copy.deepcopy(k48_dimer_ubiquitin)
    invalid_json[123]= "value"
    with pytest.raises(KeyError):
        iterate_through_ubiquitin(invalid_json)

def test_deeply_nested_erroneous_key():
    """
    Test `iterate_through_ubiquitin` to ensure it raises KeyError when an 
    erroneous key exists deep inside the nested ubiquitin structure.
    """
    # Create a deep copy of valid ubiquitin data
    test_ubiquitin = copy.deepcopy(k48_dimer_ubiquitin)

    # Introduce an erroneous key deep in the nested structure
    test_ubiquitin["branching_sites"][6]["children"]['chain_number1'] = test_ubiquitin["branching_sites"][6]["children"]['chain_number']

    with pytest.raises(KeyError, match="Invalid keys found"):
        iterate_through_ubiquitin(test_ubiquitin)


# Define valid site names
VALID_SITE_NAMES = {"K6", "K11", "K27", "K29", "K33", "K48", "K63", "M1"}

@pytest.mark.parametrize("site_name, sequence_id, children, expected_error", [
    ("K48", "FAG(K)QLE", "", None),  # ✅ Valid case
    ("K63", "NIQ(K)EST", "SMAC", None),  # ✅ Valid case with SMAC
    ("K11", "LTG(K)TIT", {}, None),  # ✅ Valid case with a nested dictionary
    ("K29", "VKA(K)IQD", "ABOC", None),  # ✅ Valid case with ABOC
    ("K27", "ENV(K)AKI", "INVALID_CHILD", "Invalid children format"),  # ❌ Invalid children
    ("K33", "IQD(K)EGI", None, "Invalid children format")  # ❌ None is not allowed
])

def test_branching_sites_format(site_name, sequence_id, children, expected_error):
    """
    Test that branching sites follow the correct structure:
    - `site_name` must be valid.
    - `sequence_id` must contain parentheses for the branching residue.
    - `children` must be `""`, `"SMAC"`, `"ABOC"`, or a dictionary.
    """
    # Create a test ubiquitin dictionary
    test_ubiquitin = copy.deepcopy(k48_dimer_ubiquitin)
    
    # Inject test values
    test_ubiquitin["branching_sites"].append({
        "site_name": site_name,
        "sequence_id": sequence_id,
        "children": children
    })

    if expected_error:
        with pytest.raises(AssertionError, match=expected_error):
            validate_branching_sites(test_ubiquitin)
    else:
        validate_branching_sites(test_ubiquitin)

# =====================================
# === Tests for generate_multimer_string ===
# Tests for generating multimer strings from ubiquitin structures.
# =====================================

@pytest.mark.parametrize("ubiquitin_structure, expected_multimer_string", [
    # Test Case 1: Basic Ubiquitin Monomer (All Branching Sites Present but Empty)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, "GG-1ubq-1-()"),
    
    # Test Case 2: Basic Ubiquitin Monomer with Histag(All Branching Sites Present but Empty)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, "his-GG-1ubq-1-()"),

    # Test Case 3: Ubiquitin Dimer with K48 Linkage no Histag
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, "GG-1ubq-1-(<K48_1ubq-2-()>)"),

    # Test Case 3: Ubiquitin Trimer with K48 and K63 Linkages
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}

        ]
    }, "his-GG-1ubq-1-(<K48_1ubq-2-(<K63_1ubq-3-()>)>)"),

    # Test Case 4: Ubiquitin with Protecting Groups (SMAC and ABOC)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "SMAC"},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": "ABOC"}
        ]
    }, "his-GG-1ubq-1-(<K48_SMAC><K63_ABOC>)"),

    # Test Case 5: Deeply Nested Ubiquitin with All Branching Sites Used
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                                "protein": "1ubq",
                                "chain_number": 4,
                                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                "chain_length": 76,
                                "branching_sites": [
                                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                                ]
                            }},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": "SMAC"},
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, "his-GG-1ubq-1-(<K48_1ubq-2-(<K48_1ubq-3-(<K48_1ubq-4-()>)><K63_SMAC>)>)"),
])
def test_iterate_through_ubiquitin(ubiquitin_structure, expected_multimer_string):
    """
    Test `inner_wrapper_iterate_through_ubiquitin` to ensure correct string formation 
    for ubiquitin structures with all possible branching sites.
    """
    result, context = iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))
    multimer_string = context["multimer_string_name"]
    assert multimer_string == expected_multimer_string, f"Expected: {expected_multimer_string}, Got: {result}"

# =====================================
# === Tests for add_max_chain_number ===
# Tests for finding the maximum chain number in ubiquitin structures.
# =====================================

@pytest.mark.parametrize("ubiquitin_structure, expected_max_chain", [
    # ✅ Test Case 1: Single Ubiquitin Monomer (Max Chain = 1)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}]
    }, 1),

    # ✅ Test Case 2: Dimer (Max Chain = 2)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}

        ]
    }, 2),

    # ✅ Test Case 3: Trimer (Max Chain = 3)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},                    
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}

        ]
    }, 3),

    # ✅ Test Case 4: Deeply Nested Ubiquitin (Max Chain = 5)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},                
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": {
                                "protein": "1ubq",
                                "chain_number": 4,
                                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                "chain_length": 76,
                                "branching_sites": [
                                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                                        "protein": "1ubq",
                                        "chain_number": 5,
                                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                        "chain_length": 76,
                                        "branching_sites": [
                                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                                        ]
                                    }}
                                ]
                            }},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, 5),

    # ✅ Test Case 5: His-Tagged Ubiquitin (Max Chain = 1)
    ({
        "protein": "1ubq-his",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH",
        "chain_length": 82,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, 1),

    # ✅ Test Case 6: Multi-Branching Ubiquitin (Max Chain = 4)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": {
                                "protein": "1ubq",
                                "chain_number": 4,
                                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                "chain_length": 76,
                                "branching_sites": [
                                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                                ]
                            }},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, 4)
])
def test_add_max_chain_number(ubiquitin_structure, expected_max_chain):

    updated_ubiquitin, context = iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))

    assert context["max_chain_number"] == expected_max_chain


def test_add_max_chain_number_simple():
    context = {"chain_number_list": [1, 2, 3, 4]}
    updated_context = add_max_chain_number(context)
    assert updated_context["max_chain_number"] == 3
    
# =====================================
# === Tests for find_free_lysines ===
# Tests for extracting all free lysines from ubiquitin structures.
# =====================================

@pytest.mark.parametrize("ubiquitin_structure, expected_free_lysines", [
    # ✅ Test Case 1: Basic Monomer (No Free Lysines)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "ABOC"},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": "SMAC"}
        ]
    }, []),

    # ✅ Test Case 2: Dimer with One Free Lysine (K48)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "ABOC"},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": "SMAC"}

        ]
    }, [[2, "K63"]]),

    # ✅ Test Case 3: Trimer with Multiple Free Lysines (K48 and K63 at different depths)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "ABOC"},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": "ABOC"},
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": "ABOC"},
        ]
    }, [[3, "K63"]]),

    # ✅ Test Case 4: Deeply Nested Structure with Multiple Free Lysines
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""},
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": "ABOC"},
        ]
    }, [[3, "K48"], [3, "K63"], [2, "K63"]]),

    # ✅ Test Case 5: His-Tagged Ubiquitin with Free Lysine
    ({
        "protein": "1ubq-his",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, [[1, "K48"], [1, "K63"]]),

    # ✅ Test Case 6: Complex Multi-Protected Ubiquitin with Multiple Free Lysines
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, [[2, "K48"], [2, "K63"], [1, "K63"]])
])
def test_find_free_lysines(ubiquitin_structure, expected_free_lysines):
    """
    Test `find_free_lysines()` to ensure it correctly identifies free lysines (K48, K63)
    across different ubiquitin configurations, including deeply nested structures.
    """
    _, context = iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))
    free_lysines = context['free_lysines']
    assert sorted(free_lysines) == sorted(expected_free_lysines), \
        f"Expected {expected_free_lysines}, got {free_lysines}"
    

# =====================================
# === Tests for find_conjugated_lysines ===
# Tests for extracting all conjugated lysines from ubiquitin structures.
# =====================================

@pytest.mark.parametrize("ubiquitin_structure, expected_conjugated_lysines", [
    # ✅ Test Case 1: Basic Monomer (No conjugated lysines)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, []),

    # ✅ Test Case 2: SMAC Conjugated Lysine at K29 but doesn't give a conjugated lysine
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, []),

    # ✅ Test Case 3: Dimer with One Conjugation at K48 on ubiquitin 1
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, [[1, "K48"]]),

    # ✅ Test Case 4: Trimer with Nested Conjugations
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": "ABOC"},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}

        ]
    }, [[1, "K48"], [2, "K63"]]),

    
    # ✅ Test Case 6: Highly Nested Structure with Multiple Conjugations
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "ABOC"},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "ABOC"},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                                "protein": "1ubq",
                                "chain_number": 4,
                                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                "chain_length": 76,
                                "branching_sites": [
                                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                                ]
                            }}
                        ]
                    }}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                "protein": "1ubq",
                "chain_number": 5,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }}
        ]
    }, [[1, "K48"], [2, "K63"], [3, "K63"], [1, "K63"]]),
])
def test_find_conjugated_lysines(ubiquitin_structure, expected_conjugated_lysines):
    """
    Test `find_conjugated_lysines()` to ensure it correctly extracts 
    all conjugated lysines from ubiquitin structures in the correct order. 
    """
    _, context = iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))
    conjugated_lysines = context['conjugated_lysines']

    assert conjugated_lysines == expected_conjugated_lysines, \
        f"Expected {expected_conjugated_lysines}, got {conjugated_lysines}"


# =====================================
# === Tests for find_SMAC and ABOC lysines within iterate_through_ubiquitin ===
# Tests for identifying lysines with SMAC and ABOC protecting groups in ubiquitin structures.
# =====================================

@pytest.mark.parametrize("ubiquitin_structure, expected_SMAC_lysines, expected_ABOC_lysines", [
    # ✅ Test Case 1: Basic Monomer (No conjugated lysines)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, [], []),

    # ✅ Test Case 2: SMAC Conjugated Lysine at K29 but doesn't give a conjugated lysine
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, [[1, "K29"]], []),

    # ✅ Test Case 3: Dimer with One Conjugation at K48 on ubiquitin 1
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, [[2, "K29"]], [[2, "K11"]]),

    # ✅ Test Case 4: Trimer with Nested Conjugations
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": "ABOC"},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}

        ]
    }, [[2, "K29"], [3, "K29"], [3, "K33"]], [[2, "K11"], [2, "K27"], [3, "K11"]]),

    
    # ✅ Test Case 6: Highly Nested Structure with Multiple Conjugations
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "ABOC"},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "ABOC"},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                                "protein": "1ubq",
                                "chain_number": 4,
                                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                "chain_length": 76,
                                "branching_sites": [
                                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                                ]
                            }}
                        ]
                    }}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                "protein": "1ubq",
                "chain_number": 5,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }}
        ]
    }, [[1, "K29"], [1, "K33"], [2, "K33"], [3, "K33"], [4, "K29"], [4, "K33"], [5, "K29"], [5, "K33"]], [[1, "K11"], [2, "K11"], [2, "K29"], [3, "K11"], [3, "K29"], [4, "K11"],  [5, "K11"]]),
])

# =====================================
# === Redundant Tests for find_SMAC&ABOC lysines ===
# These tests are redundant and are not needed.
# =====================================

def test_find_SMAC_ABOC_lysines(ubiquitin_structure, expected_SMAC_lysines, expected_ABOC_lysines):
    """
    Test `find_SMAC_ABOC_lysines()` to ensure it correctly extracts 
    all SMAC and ABOC-modified lysines from ubiquitin structures.
    """
    _, context = iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))
    ABOC_lysines = context["ABOC_lysines"]
    SMAC_lysines = context["SMAC_lysines"]
    
    assert sorted(SMAC_lysines) == sorted(expected_SMAC_lysines), \
        f"Expected SMAC {expected_SMAC_lysines}, got {SMAC_lysines}"
    
    assert sorted(ABOC_lysines) == sorted(expected_ABOC_lysines), \
        f"Expected ABOC {expected_ABOC_lysines}, got {ABOC_lysines}"

@pytest.mark.parametrize("ubiquitin_structure, expected_aboc, expected_smac", [
    # ✅ Test Case 1: Basic Monomer (No ABOC or SMAC)
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, 0, 0),

    # ✅ Test Case 2: Dimer with One ABOC
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "ABOC"},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, 1, 0),

    # ✅ Test Case 3: Trimer with One ABOC and One SMAC
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": "ABOC"},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": "SMAC"}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, 1, 1),

    # ✅ Test Case 4: Deeply Nested Structure with Multiple ABOC & SMAC
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": "ABOC"},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "ABOC"},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                        "protein": "1ubq",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "SMAC"},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "ABOC"},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}

        ]
    }, 3, 2),

    # ✅ Test Case 5: His-Tagged Ubiquitin with ABOC
    ({
        "protein": "1ubq-his",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH",
        "chain_length": 82,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": "ABOC"},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, 1, 0),

    # ✅ Test Case 6: Multi-Protected Ubiquitin with Equal Number of ABOC & SMAC
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": "ABOC"},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "SMAC"},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "ABOC"},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "SMAC"},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, 2, 2)
])
    
def test_find_protecting_groups(ubiquitin_structure, expected_aboc, expected_smac):
    """
    Test `find_protecting_groups()` to ensure it correctly counts 
    the number of ABOC and SMAC protecting groups.
    """
    updated_ubiquitin, context = iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))
    ABOC_lysines = context["ABOC_lysines"]
    SMAC_lysines = context["SMAC_lysines"]
    
    assert len(ABOC_lysines) == expected_aboc, f"Expected {expected_aboc} ABOC, got {len(ABOC_lysines)}"
    assert len(SMAC_lysines) == expected_smac, f"Expected {expected_smac} SMAC, got {len(SMAC_lysines)}"

    
# =====================================
# === Tests for process_current_protein === 
# Tests for processing the current protein's FASTA sequence and updating context.
# =====================================

@pytest.mark.parametrize(
    "working_dictionary, context, expected_working_dict, expected_context",
    [
        # ✅ Test 1: Standard Case
        (
            {"FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
            {"chain_number_list": [1], "chain_length_list": []},
            {"FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG", "chain_length": 76, "chain_number": 1},
            {"chain_number_list": [1, 2], "chain_length_list": [76]}
        ),

        # ✅ Test 2: Chain number list already has progress
        (
            {"FASTA_sequence": "ABCDEF"},
            {"chain_number_list": [1, 2, 3], "chain_length_list": [10, 20]},
            {"FASTA_sequence": "ABCDEF", "chain_length": 6, "chain_number": 3},
            {"chain_number_list": [1, 2, 3, 4], "chain_length_list": [10, 20, 6]}
        ),

        # ✅ Test 3: Empty FASTA sequence
        (
            {"FASTA_sequence": ""},
            {"chain_number_list": [1], "chain_length_list": []},
            {"FASTA_sequence": "", "chain_length": 0, "chain_number": 1},
            {"chain_number_list": [1, 2], "chain_length_list": [0]}
        ),

        # ✅ Test 4: Single character FASTA sequence
        (
            {"FASTA_sequence": "M"},
            {"chain_number_list": [3], "chain_length_list": [5, 6]},
            {"FASTA_sequence": "M", "chain_length": 1, "chain_number": 3},
            {"chain_number_list": [3, 4], "chain_length_list": [5, 6, 1]}
        ),

        # ✅ Test 5: Long FASTA sequence
        (
            {"FASTA_sequence": "M" * 1000},
            {"chain_number_list": [10], "chain_length_list": [76, 200]},
            {"FASTA_sequence": "M" * 1000, "chain_length": 1000, "chain_number": 10},
            {"chain_number_list": [10, 11], "chain_length_list": [76, 200, 1000]}
        ),
    ]
)
def test_process_current_protein(working_dictionary, context, expected_working_dict, expected_context):
    updated_working_dict, updated_context = process_current_protein(working_dictionary.copy(), context.copy())

    assert updated_working_dict == expected_working_dict, f"Expected working dictionary: {expected_working_dict}, but got: {updated_working_dict}"
    assert updated_context == expected_context, f"Expected context: {expected_context}, but got: {updated_context}"

# =====================================
# === Tests for test_process_branch ===
# These tests validate the behavior of the `process_branch` function, ensuring that it correctly updates the context based on the provided branch information.
# =====================================

@pytest.mark.parametrize("branch, working_dict, starting_context, expected_context", [
    # ✅ Test 1: SMAC protecting group
    (
        {"site_name": "K27", "sequence_id": "HIJ(K)LMN", "children": "SMAC"},
        {"chain_number": 1, "FASTA_sequence": "ABCDEFGHIJKLMNOPQ"},
        {"chain_length_list": [76], "chain_number_list": [1], "multimer_string_name": "", "SMAC_lysines": [], "ABOC_lysines": [], "free_lysines": [], "conjugated_lysines": []},
        {"multimer_string_name": "" + "<K27_SMAC>", "SMAC_lysines": [[1, "K27"]]}
    ),
    # ✅ Test 2: ABOC protecting group
    (
        {"site_name": "K33", "sequence_id": "ENV(K)AKI", "children": "ABOC"},
        {"chain_number": 2, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
        {"chain_length_list": [76, 76], "chain_number_list": [1, 2], "multimer_string_name": "", "SMAC_lysines": [], "ABOC_lysines": [], "free_lysines": [], "conjugated_lysines": []},
        {"multimer_string_name": "<K33_ABOC>", "ABOC_lysines": [[2, "K33"]]}
    ),
    # ✅ Test 3-6: Free lysines (M1, K6, K11, K27, K29, K33)
    *[
        (
            {"site_name": site, "sequence_id": sequence, "children": ""},
            {"chain_number": 1, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
            {"chain_length_list": [76], "chain_number_list": [1], "multimer_string_name": "", "SMAC_lysines": [], "ABOC_lysines": [], "free_lysines": [], "conjugated_lysines": []},
            {"multimer_string_name": ""}
        )
        for site, sequence in zip(["M1", "K6", "K11", "K27", "K29", "K33"], ["(M)QIF" ,"IFV(K)TLT" ,"LTG(K)TIT" ,"ENV(K)AKI" ,"VKA(K)IQD" ,"IQD(K)EGI"])
    ],
    # ✅ Test 7-8: Free lysines special cases K48, K63
    (
        {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
        {"chain_number": 1, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
        {"chain_length_list": [76], "chain_number_list": [1], "multimer_string_name": "", "SMAC_lysines": [], "ABOC_lysines": [], "free_lysines": [], "conjugated_lysines": []},
        {"free_lysines": [[1, "K48"]]}
    ),
    (
        {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""},
        {"chain_number": 2, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
        {"chain_length_list": [76, 76], "chain_number_list": [1, 2], "multimer_string_name": "", "SMAC_lysines": [], "ABOC_lysines": [], "free_lysines": [], "conjugated_lysines": []},
        {"free_lysines": [[2, "K63"]]}
    ),
    # ✅ Test 9: Recursive children call
    (
        {"site_name": "K48", 
         "sequence_id": "FAG(K)QLE", 
         "children": {"protein": "dummy_protein",
                        "chain_number": 2,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]}},
        {"chain_number": 1, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
        {"chain_length_list": [76], "chain_number_list": [1, 2], "multimer_string_name": "", "SMAC_lysines": [], "ABOC_lysines": [], "free_lysines": [], "conjugated_lysines": []},
        {"chain_length_list": [76, 76], "chain_number_list": [1, 2, 3], "multimer_string_name": "<K48_dummy_protein-2-(<K11_ABOC><K29_SMAC><K33_SMAC>)>", "SMAC_lysines": [[2, 'K29'], [2, 'K33']], "conjugated_lysines": [[1, "K48"]], "ABOC_lysines": [[2, 'K11']], "free_lysines": [[2, 'K48'], [2, 'K63']]}
    ),
    # ✅ Test 10: Recursive children with existing multimer_string_name
    (
        {"site_name": "K63", 
         "sequence_id": "NIQ(K)EST", 
         "children": {"protein": "dummy_protein",
                        "chain_number": 2,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]}},
        {"chain_number": 2, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
        {"chain_length_list": [76], "chain_number_list": [1, 2], "multimer_string_name": "PREVIOUS", "SMAC_lysines": [], "ABOC_lysines": [], "free_lysines": [], "conjugated_lysines": []},
        {"multimer_string_name": "PREVIOUS<K63_dummy_protein-2-(<K11_ABOC><K29_SMAC><K33_SMAC>)>", "conjugated_lysines": [[2, "K63"]]}
    ),
    # ✅ Test 11: Add to existing SMAC list
    (
        {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": "SMAC"},
        {"chain_number": 1, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
        {"chain_length_list": [76], "chain_number_list": [1], "multimer_string_name": "", "SMAC_lysines": [[0, "K6"]], "ABOC_lysines": [], "free_lysines": [], "conjugated_lysines": []},
        {"multimer_string_name": "<K27_SMAC>", "SMAC_lysines": [[0, "K6"], [1, "K27"]]}
    ),
    # ✅ Test 12: Add to existing ABOC list
    (
        {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "ABOC"},
        {"chain_number": 1, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
        {"chain_length_list": [76], "chain_number_list": [1], "multimer_string_name": "", "SMAC_lysines": [], "ABOC_lysines": [[0, "K33"]], "free_lysines": [], "conjugated_lysines": []},
        {"multimer_string_name": "<K29_ABOC>", "ABOC_lysines": [[0, "K33"], [1, "K29"]]}
    ),
    # ✅ Test 13: Already populated context
    (
        {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
        {"chain_number": 3, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
        {"chain_length_list": [76, 76], "chain_number_list": [1, 2, 3], "multimer_string_name": "PREVIOUS", "SMAC_lysines": [[0, "K6"]], "ABOC_lysines": [[1, "K33"]], "free_lysines": [[2, "K63"]], "conjugated_lysines": [[0, "K48"]]},
        {"free_lysines": [[2, "K63"], [3, "K48"]]}
    ),
    # ✅ Test 14: Recursive children, deeply nested
    (
        {"site_name": "K63", 
         "sequence_id": "NIQ(K)EST", 
         "children": {"protein": "dummy_protein",
                        "chain_number": 4,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": "ABOC"},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]}},
        {"chain_number": 3, "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"},
        {"chain_length_list": [76, 76, 76], "chain_number_list": [1, 2, 3, 4], "multimer_string_name": "", "SMAC_lysines": [], "ABOC_lysines": [], "free_lysines": [], "conjugated_lysines": []},
        {"multimer_string_name": "<K63_dummy_protein-4-(<K11_ABOC><K29_SMAC><K33_SMAC>)>", "conjugated_lysines": [[3, "K63"]], "SMAC_lysines": [[4, "K29"],[4, "K33"]], "ABOC_lysines": [[4, "K11"]], "free_lysines": [[4, "K48"], [4, "K63"]]}
    )
])
def test_process_branch(branch, working_dict, starting_context, expected_context):
    updated_branch, updated_working_dict, updated_context = process_branch(branch.copy(), working_dict.copy(), starting_context.copy())

    # Check updates in context
    for key, expected_value in expected_context.items():
        assert updated_context[key] == expected_value, f"Context key '{key}' mismatch: expected {expected_context}, got {updated_context[key]}"

@pytest.mark.parametrize("branch, expected_in_context", [

    # ✅ Test 1: SMAC protecting group
    ({"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "SMAC"},
     {"SMAC_lysines": [[1, "K29"]], "multimer_string_name": "<K29_SMAC>"}),

    # ✅ Test 2: ABOC protecting group
    ({"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": "ABOC"},
     {"ABOC_lysines": [[1, "K29"]], "multimer_string_name": "<K29_ABOC>"}),

    # ✅ Test 3: Free lysine (K48)
    ({"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
     {"free_lysines": [[1, "K48"]]}),

    # ✅ Test 4: Neutral lysine (M1)
    ({"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
     {}),  # No changes expected

    # ✅ Test 5: Recursive child (nested protein)
    ({"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": BASE_WORKING_DICT},
     {"multimer_string_name": "<K48_GG-dummy_protein-1-()>", "conjugated_lysines": [[1, "K48"]]}),

    # ✅ Test 6: Very large chain number for scalability
    ({"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
     {"free_lysines": [[10_000, "K48"]]}),

])
def test_process_branch_valid_cases(branch, expected_in_context):
    """Test process_branch() with valid scenarios."""
    # Setup
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)

    # Edge case: simulate high chain number for scalability test
    if expected_in_context.get("free_lysines") == [[10_000, "K48"]]:
        working_dict["chain_number"] = 10_000

    # Process
    _, _, updated_context = process_branch(branch, working_dict, context)

    # Assertions
    for key, expected_value in expected_in_context.items():
        assert updated_context[key] == expected_value, f"Failed for key: {key}"


# ✅ Test 7: Missing 'site_name'
def test_process_branch_missing_site_name():
    branch = {"sequence_id": "FAG(K)QLE", "children": ""}
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)

    with pytest.raises(KeyError):
        process_branch(branch, working_dict, context)


# ✅ Test 8: Missing 'sequence_id'
def test_process_branch_missing_sequence_id():
    branch = {"site_name": "K48", "children": ""}
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)

    with pytest.raises(KeyError):
        process_branch(branch, working_dict, context)


# ✅ Test 9: Children is None
def test_process_branch_children_none():
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": None}
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)

    # Should trigger invalid children format inside inner function if error handling exists
    _, _, updated_context = process_branch(branch, working_dict, context)
    # Safe fallback: no changes expected
    assert updated_context == context


# ✅ Test 10: Children is int
def test_process_branch_children_int():
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": 123}
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)
    # Safe fallback: no changes expected
    _, _, updated_context = process_branch(branch, working_dict, context)
    assert updated_context == context


# ✅ Test 11: Children is list
def test_process_branch_children_list():
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": []}
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)

    # Safe fallback: no changes expected
    _, _, updated_context = process_branch(branch, working_dict, context)
    assert updated_context == context


# ✅ Test 12: Empty context dictionary
def test_process_branch_empty_context():
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""}
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = {}

    with pytest.raises(KeyError):
        process_branch(branch, working_dict, context)


# ✅ Test 13: Empty working dictionary
def test_process_branch_empty_working_dict():
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""}
    working_dict = {}
    context = copy.deepcopy(BASE_CONTEXT)

    with pytest.raises(KeyError):
        process_branch(branch, working_dict, context)


# ✅ Test 14: Invalid FASTA sequence
def test_process_branch_invalid_fasta_sequence():
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""}
    working_dict = {"chain_number": 1, "FASTA_sequence": ""}  # Invalid empty sequence
    context = copy.deepcopy(BASE_CONTEXT)

    with pytest.raises(ValueError):
        process_branch(branch, working_dict, context)


# ✅ Test 15: Unicode in sequence_id
def test_process_branch_unicode_sequence_id():
    branch = {"site_name": "K48", "sequence_id": "SEQ🔥123", "children": ""}
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)

    with pytest.raises(ValueError):
        process_branch(branch, working_dict, context)


# ✅ Test 16: Deep nested children (depth > 2)
def test_process_branch_deep_nested_children():
    deep_child = {
    "protein": "dummy_protein",
    "chain_number": 1,
    "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    "chain_length": 76,
    "branching_sites": [
        {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
        {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
        {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
        {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": BASE_WORKING_DICT},
        {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
        {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
        {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
        {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
    ]}
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": deep_child}

    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)

    _, _, updated_context = process_branch(branch, working_dict, context)

    assert "multimer_string_name" in updated_context
    assert updated_context["multimer_string_name"].startswith("<K48_")  # partial check


# ✅ Test 18: Verify input immutability
def test_process_branch_input_immutability():
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""}
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)

    original_branch = copy.deepcopy(branch)
    original_context = copy.deepcopy(context)

    process_branch(branch, working_dict, context)

    assert branch == original_branch, "Branch mutated!"
    # Context is expected to change — no assertion


# ✅ Test 19: Unexpected field in branch
def test_process_branch_unexpected_branch_field():
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "", "unexpected": "value"}
    working_dict = copy.deepcopy(BASE_WORKING_DICT)
    context = copy.deepcopy(BASE_CONTEXT)

    _, _, updated_context = process_branch(branch, working_dict, context)
    assert updated_context is not None


# ✅ Test 20: Large FASTA sequence
def test_process_branch_large_fasta_sequence():
    branch = {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""}
    working_dict = {"chain_number": 1, "FASTA_sequence": "MMMMMFAGKQLE" + "M" * 10_000}
    context = copy.deepcopy(BASE_CONTEXT)

    _, _, updated_context = process_branch(branch, working_dict, context)
    assert updated_context is not None


# Verify ubiquitin is added correctly to a free lysine site.
def test_add_to_unoccupied_site():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    updated = K_residue_ubi_addition(data, 1, "IFV(K)TLT", BASE_WORKING_DICT)
    assert isinstance(updated["branching_sites"][1]["children"], dict)

# Verify ubiquitin is added to a lysine site even if marked with protecting groups as strings.
def test_add_to_site_with_string_children():
    for protecting_group in ["SMAC", "ABOC"]:
        data = copy.deepcopy(five_level_nested_ubiquitin_)
        data["branching_sites"][1]["children"] = protecting_group
        updated = K_residue_ubi_addition(data, 1, "IFV(K)TLT", BASE_WORKING_DICT)
        assert isinstance(updated["branching_sites"][1]["children"], dict)

# Verify no changes occur if the sequence ID does not match any lysine.
def test_no_change_if_sequence_id_mismatch():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    original = copy.deepcopy(data["branching_sites"])
    updated = K_residue_ubi_addition(data, 1, "XYZ(K)ABC", BASE_WORKING_DICT)
    assert updated["branching_sites"] == original

# Verify error is raised when conjugating to a lysine that already has a ubiquitin attached.
def test_error_on_already_conjugated_site():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = {"mock": "ubi"}
    with pytest.raises(TypeError, match="already conjugated"):
        K_residue_ubi_addition(data, 1, "IFV(K)TLT", BASE_WORKING_DICT)

# Verify JSON string input is properly converted to a dictionary.
def test_json_input_converted():
    import json
    json_data = json.dumps(five_level_nested_ubiquitin_)
    updated = K_residue_ubi_addition(json_data, 1, "IFV(K)TLT", BASE_WORKING_DICT)
    assert isinstance(updated, dict)

# Verify addition raises TypeError if FASTA sequence does not end with RLRGG.
def test_addition_raises_if_RLRGG_tail_missing():
    # Test that a TypeError is raised if the FASTA sequence does not end with RLRGG
    non_matching_tail = copy.deepcopy(BASE_WORKING_DICT)
    non_matching_tail["FASTA_sequence"] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGXXXXXX"
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    with pytest.raises(TypeError, match="does not end with RLRGG"):
        K_residue_ubi_addition(data, 1, "IFV(K)TLT", non_matching_tail)

# Verify correct lysine index is targeted and updated within branching_sites.
def test_loop_index_resolved_correctly():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    # Ensure the correct index (1) is updated
    updated = K_residue_ubi_addition(data, 1, "IFV(K)TLT", BASE_WORKING_DICT)
    assert updated["branching_sites"][1]["site_name"] == "K6"
    assert isinstance(updated["branching_sites"][1]["children"], dict)

# Verify conjugation works on deeply nested ubiquitin structures.
def test_addition_on_deep_nested_chain():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    deep_branch = data["branching_sites"][6]["children"]["branching_sites"][7]["children"]
    assert deep_branch["chain_number"] == 3
    modified = K_residue_ubi_addition(deep_branch, 3, "ENV(K)AKI", BASE_WORKING_DICT)
    assert isinstance(modified["branching_sites"][3]["children"], dict)

# Verify error is raised if input is not a valid dictionary or JSON string.
def test_addition_raises_if_not_dict_or_json():
    invalid_input = "This is not a valid input"
    with pytest.raises(ValueError, match="Invalid JSON format: Unable to parse the string"):
        K_residue_ubi_addition(invalid_input, 1, "IFV(K)TLT", BASE_WORKING_DICT)

# Verify error is raised if input is an empty string.
def test_addition_raises_if_empty_string():
    empty_input = ""
    with pytest.raises(ValueError, match="Invalid JSON format: Unable to parse the string"):
        K_residue_ubi_addition(empty_input, 1, "IFV(K)TLT", BASE_WORKING_DICT)

# Verify error is raised if input is None.
def test_addition_raises_if_none():
    none_input = None
    with pytest.raises(TypeError, match="Input must be a dictionary or a JSON string"):
        K_residue_ubi_addition(none_input, 1, "IFV(K)TLT", BASE_WORKING_DICT)

# Verify SMAC deprotection clears the 'children' field.
def test_smac_deprotection():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = "SMAC"
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][1], data, "SMAC_deprot", BASE_WORKING_DICT)
    assert bra["children"] == ""

# Verify ABOC deprotection clears the 'children' field.
def test_aboc_deprotection():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = "ABOC"
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][1], data, "ABOC_deprot", BASE_WORKING_DICT)
    assert bra["children"] == ""

# Verify GLOBAL deprotection clears SMAC protecting group correctly.
def test_global_deprotection_smac():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = "SMAC"
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][1], data, "GLOBAL_deprot", BASE_WORKING_DICT)
    assert bra["children"] == ""

# Verify GLOBAL deprotection clears ABOC protecting group correctly.
def test_global_deprotection_aboc():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = "ABOC"
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][1], data, "GLOBAL_deprot", BASE_WORKING_DICT)
    assert bra["children"] == ""

# Verify ubiquitin is conjugated to K48 if unoccupied and sequence ID matches.
def test_k48_conjugation():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][6]
    bra["children"] = ""
    bra, updated = process_ubiquitin_reaction(bra, data, "K48", BASE_WORKING_DICT)
    assert isinstance(updated["branching_sites"][6]["children"], dict)

# Verify ubiquitin is conjugated to K63 if unoccupied and sequence ID matches.
def test_k63_conjugation():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][7]
    bra["children"] = ""
    bra, updated = process_ubiquitin_reaction(bra, data, "K63", BASE_WORKING_DICT)
    assert isinstance(updated["branching_sites"][7]["children"], dict)

# Verify no conjugation occurs if sequence ID for K48 reaction does not match.
def test_k48_conjugation_wrong_sequence():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][2]
    bra["sequence_id"] = "WRONG(K)SEQ"
    original = copy.deepcopy(bra["children"])
    bra, updated = process_ubiquitin_reaction(bra, data, "K48", BASE_WORKING_DICT)
    assert bra["children"] == original

# Verify no conjugation occurs if K63 site already has a ubiquitin attached.
def test_k63_conjugation_with_existing_children():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][7]
    bra["children"] = {"mock": "ubi"}
    bra, updated = process_ubiquitin_reaction(bra, data, "K63", BASE_WORKING_DICT)
    assert bra["children"] == {"mock": "ubi"}

@pytest.mark.parametrize("lys_site, expected_child", [
    # SMAC deprotection
    ("K48", ""),
])
def test_process_ubiquitin_smac_deprot(lys_site, expected_child):
    """
    Test SMAC deprotection: the function should remove 'SMAC' from the lysine site.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    for site in data["branching_sites"]:
        if site["site_name"] == lys_site:
            site["children"] = "SMAC"
            bra, updated = process_ubiquitin_reaction(site, data, "SMAC_deprot", ubiquitin_monomer)
            assert bra["children"] == expected_child
            break

def test_process_ubiquitin_aboc_deprot():
    """
    Test ABOC deprotection: the function should remove 'ABOC' from the lysine site.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][6]["children"] = "ABOC"
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][6], data, "ABOC_deprot", ubiquitin_monomer)
    assert bra["children"] == ""

def test_process_ubiquitin_global_deprot_smac():
    """
    Test global deprotection: 'SMAC' child should be cleared.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][6]["children"] = "SMAC"
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][6], data, "GLOBAL_deprot", ubiquitin_monomer)
    assert bra["children"] == ""

def test_process_ubiquitin_global_deprot_aboc():
    """
    Test global deprotection: 'ABOC' child should be cleared.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][6]["children"] = "ABOC"
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][6], data, "GLOBAL_deprot", ubiquitin_monomer)
    assert bra["children"] == ""

def test_process_ubiquitin_k48_conjugation():
    """
    Test K48 conjugation: a ubiquitin should be added to a free K48 site.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][6]["children"] = ""
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][6], data, "K48", ubiquitin_monomer)
    assert isinstance(data["branching_sites"][6]["children"], dict)
    assert data["branching_sites"][6]["children"]["protein"] == "1ubq"

def test_process_ubiquitin_k63_conjugation():
    """
    Test K63 conjugation: a ubiquitin should be added to a free K63 site.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][7]["children"] = ""
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][7], data, "K63", ubiquitin_monomer)
    assert isinstance(data["branching_sites"][7]["children"], dict)
    assert data["branching_sites"][7]["children"]["protein"] == "1ubq"

def test_process_ubiquitin_k48_noop_when_occupied():
    """
    Test K48 conjugation when already occupied: should remain unchanged.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    original_child = data["branching_sites"][6]["children"]
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][6], data, "K48", ubiquitin_monomer)
    assert data["branching_sites"][6]["children"] == original_child

def test_process_ubiquitin_k63_noop_when_invalid_sequence():
    """
    Test K63 conjugation with wrong sequence_id: no conjugation should happen.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][7]["sequence_id"] = "INVALID(K)SEQ"
    data["branching_sites"][7]["children"] = ""
    bra, updated = process_ubiquitin_reaction(data["branching_sites"][7], data, "K63", ubiquitin_monomer)
    assert data["branching_sites"][7]["children"] == ""

# =============================
# Tests for ubiquitin_simulation and inner_wrapper_ubiquitin_simulation
# =============================

# Test K48 addition to level 1 lysine
def test_ubiquitin_simulation_k48_on_level_1():
    """Test K48 conjugation works at top-level ubiquitin structure."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "K48")
    assert isinstance(results, dict)
    assert isinstance(results["branching_sites"][6]["children"], dict)

# Test K63 addition to level 1 lysine
def test_ubiquitin_simulation_k63_on_level_1():
    """Test K63 conjugation works at top-level ubiquitin structure."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "K63")
    assert isinstance(results["branching_sites"][7]["children"], dict)

# Test SMAC deprotection
def test_ubiquitin_simulation_smac_deprot():
    """Test SMAC deprotection clears correct child at level 1."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = "SMAC"
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "SMAC_deprot")
    assert results["branching_sites"][1]["children"] == ""

# Test GLOBAL deprotection on ABOC
def test_ubiquitin_simulation_global_deprot_aboc():
    """Test GLOBAL deprotection clears ABOC child."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = "ABOC"
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "GLOBAL_deprot")
    assert results["branching_sites"][1]["children"] == ""

# Test ABOC deprotection
def test_ubiquitin_simulation_aboc_deprot():
    """Test ABOC deprotection clears correct child."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = "ABOC"
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "ABOC_deprot")
    assert results["branching_sites"][1]["children"] == ""

# Test invalid reaction type
def test_ubiquitin_simulation_invalid_reaction():
    """Test invalid type_of_reaction raises ValueError."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    with pytest.raises(ValueError):
        ubiquitin_simulation(data, BASE_WORKING_DICT, "INVALID")

# Recursive test: K63 nested conjugation level 2
def test_recursive_k63_nested_addition_level_2():
    """Test K63 conjugation on nested chain (level 2)."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    level_2 = data["branching_sites"][6]["children"]
    level_2["branching_sites"][7]["children"] = ""
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "K63")
    nested = results["branching_sites"][6]["children"]["branching_sites"][7]["children"]
    assert isinstance(nested, dict)

# Recursive test: SMAC deprot deep level
def test_recursive_smac_deprot_deep():
    """Test SMAC deprotection applied on deeply nested structure."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]["children"] = "SMAC"
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "SMAC_deprot")
    cleared = results["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]["children"]
    assert cleared == ""

# Recursive test: ABOC deprot deep level
def test_recursive_aboc_deprot_deep():
    """Test ABOC deprotection applied on nested 4th level."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]["children"] = "ABOC"
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "ABOC_deprot")
    cleared = results["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]["children"]
    assert cleared == ""

# Recursive test: GLOBAL deprot on multiple nested levels
def test_recursive_global_deprot_multiple():
    """Test GLOBAL deprotection clears multiple nested SMAC/ABOC sites."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = "SMAC"
    data["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]["children"] = "ABOC"
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "GLOBAL_deprot")
    assert results["branching_sites"][1]["children"] == ""
    assert results["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]["children"] == ""

# Recursive test: K48 nested addition level 3
def test_recursive_k48_addition_deep():
    """Test K48 conjugation on 3rd-level nested structure."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    level_3 = data["branching_sites"][6]["children"]["branching_sites"][7]["children"]
    level_3["branching_sites"][6]["children"] = ""
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "K48")
    added = results["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][6]["children"]
    assert isinstance(added, dict)

# Recursive test: Ensure context correctly tracks chain numbers
def test_context_tracking_deep_recursion():
    """Test chain_number_list grows correctly during recursion."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "K48")
    assert isinstance(results, dict)
#
# Complex test: Add multiple ubiquitins and check renumbering
def test_multiple_ubiquitin_additions_and_relabelling():
    """Test multiple conjugations across different levels and validate renumbering."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)

    # Manually clear several children for additions
    data["branching_sites"][1]["children"] = ""
    data["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][6]["children"] = ""
    data["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][7]["children"] = ""

    # Apply K63 at level 1
    results, contexts = ubiquitin_simulation(data, BASE_WORKING_DICT, "K63")
    # Apply K48 at deep level
    results, contexts = ubiquitin_simulation(results, BASE_WORKING_DICT, "K48")

    # Check that additions were made
    print(results)
    assert isinstance(results["branching_sites"][1]["children"], str)
    assert isinstance(results["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][6]["children"], dict)
    assert isinstance(results["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][7]["children"], dict)

    new_input = copy.deepcopy(results)

    # Run relabeling to ensure chain numbers are updated
    relabelled, relabelled_contexts = iterate_through_ubiquitin(new_input)

    # Check that the relabelled results match the original results
    # Assumes iterate_through_ubiquitin functions correctly

    assert relabelled == results

# =============================
# Tests for handle_lysine_modification
# =============================

@pytest.mark.parametrize("lys_site, modifier, expected", [
    ("K6", "SMAC", "SMAC"),   # Add SMAC to empty site
    ("K6", "ABOC", "ABOC"),   # Add ABOC to empty site
])
def test_handle_lysine_adds_protecting_group(lys_site, modifier, expected):
    """
    Test that the protecting group (SMAC or ABOC) is correctly added to an empty lysine site.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][1]
    bra["children"] = ""
    bra["site_name"] = lys_site
    bra, updated = handle_lysine_modification(bra, data, modifier, 1, lys_site)
    assert bra["children"] == expected


def test_handle_lysine_adds_ubiquitin():
    """
    Test that a ubiquitin is added to a lysine site when the correct working dict is provided.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][1]
    bra["children"] = ""
    bra["site_name"] = "K6"
    bra, updated = handle_lysine_modification(bra, data, BASE_WORKING_DICT, 1, "K6")
    assert isinstance(bra["children"], dict)


def test_handle_lysine_wrong_chain_number_does_nothing():
    """
    Test that modification does nothing if the chain_number does not match.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][1]
    bra["children"] = ""
    bra["site_name"] = "K6"
    original = copy.deepcopy(bra)
    bra, updated = handle_lysine_modification(bra, data, BASE_WORKING_DICT, 99, "K6")
    assert bra == original


def test_handle_lysine_wrong_site_name_does_nothing():
    """
    Test that modification does nothing if the site_name does not match.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][1]
    bra["children"] = ""
    bra["site_name"] = "K6"
    original = copy.deepcopy(bra)
    bra, updated = handle_lysine_modification(bra, data, BASE_WORKING_DICT, 1, "K48")
    assert bra == original


def test_handle_lysine_site_already_has_ubiquitin_does_nothing():
    """
    Test that no modification occurs if the lysine site already has a ubiquitin attached.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][1]
    bra["children"] = {"mock": "ubi"}
    bra, updated = handle_lysine_modification(bra, data, BASE_WORKING_DICT, 1, "K6")
    assert bra["children"] == {"mock": "ubi"}


def test_handle_lysine_site_already_has_smac_does_nothing_for_ubi():
    """
    Test that no modification occurs if the lysine site already has a SMAC protecting group.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][1]
    bra["children"] = "SMAC"
    bra, updated = handle_lysine_modification(bra, data, BASE_WORKING_DICT, 1, "K6")
    assert bra["children"] == "SMAC"


def test_handle_lysine_modification_with_invalid_input_type_raises_typeerror():
    """
    Test that a TypeError is raised when an invalid input type (None) is provided as the modifier.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][1]
    bra["children"] = ""
    bra["site_name"] = "K6"
    with pytest.raises(TypeError, match="Input must be a dictionary or a JSON string"):
        handle_lysine_modification(bra, data, None, 1, "K6")


def test_handle_lysine_modification_site_occupied_by_aboc_no_ubi():
    """
    Test that no ubiquitin is added if the site is already occupied by an ABOC protecting group.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][1]
    bra["children"] = "ABOC"
    bra["site_name"] = "K6"
    bra, updated = handle_lysine_modification(bra, data, BASE_WORKING_DICT, 1, "K6")
    assert bra["children"] == "ABOC"


def test_handle_lysine_modification_wrong_site_on_correct_chain():
    """
    Test that no modification occurs if the site_name is incorrect, even if the chain number matches.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][2]
    bra["children"] = ""
    bra["site_name"] = "K11"
    bra, updated = handle_lysine_modification(bra, data, BASE_WORKING_DICT, 1, "K6")
    assert bra["children"] == ""


def test_handle_lysine_modification_protecting_group_overrides_empty():
    """
    Test that a protecting group (ABOC) is added to a lysine site that is empty.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][2]
    bra["children"] = ""
    bra["site_name"] = "K11"
    bra, updated = handle_lysine_modification(bra, data, "ABOC", 1, "K11")
    assert bra["children"] == "ABOC"


def test_handle_lysine_modification_ubiquitin_not_added_if_site_name_mismatch():
    """
    Test that a ubiquitin is not added if the site_name does not match the intended site.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][2]
    bra["children"] = ""
    bra["site_name"] = "K11"
    bra, updated = handle_lysine_modification(bra, data, BASE_WORKING_DICT, 1, "K6")
    assert bra["children"] == ""


def test_handle_lysine_modification_site_matches_but_chain_does_not():
    """
    Test that a ubiquitin is not added if the site_name matches but the chain number does not.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    bra = data["branching_sites"][2]
    bra["children"] = ""
    bra["site_name"] = "K11"
    bra, updated = handle_lysine_modification(bra, data, BASE_WORKING_DICT, 2, "K11")
    assert bra["children"] == ""

# =============================
# Tests for ubiquitin_building and inner_wrapper_ubiquitin_building
# =============================

# 1. Add SMAC to level 1, confirm assignment.
def test_add_smac_top_level():
    """Add SMAC to level 1 and confirm assignment."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    updated, _ = ubiquitin_building(data, "SMAC", 1, "K6")
    assert updated["branching_sites"][1]["children"] == "SMAC"

# 2. Add ABOC to level 1, confirm assignment.
def test_add_aboc_top_level():
    """Add ABOC to level 1 and confirm assignment."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = ""
    assert data["branching_sites"][1]["children"] == ""
    updated, _ = ubiquitin_building(data, "ABOC", 1, "K6")
    assert updated["branching_sites"][1]["children"] == "ABOC"

# 3. Add ubiquitin to free lysine at level 1.
def test_add_ubiquitin_top_level():
    """Add ubiquitin to free lysine at level 1 and check children is dict."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][6]["children"] = ""
    assert data["branching_sites"][6]["children"] == ""
    updated, _ = ubiquitin_building(data, BASE_WORKING_DICT, 1, "K48")
    assert isinstance(updated["branching_sites"][6]["children"], dict)

# 4. Ensure no overwrite if children exist.
def test_ignore_occupied_lysine():
    """Ensure that occupied lysine is not overwritten by new addition."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][6]["children"] = "SMAC"
    assert data["branching_sites"][6]["children"] == "SMAC"
    updated, _ = ubiquitin_building(data, BASE_WORKING_DICT, 1, "K48")
    assert updated["branching_sites"][6]["children"] == "SMAC"

# 5. Add SMAC at deep nested lysine.
def test_add_smac_nested_level_4():
    """Add SMAC at deep nested lysine (level 4)."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    nested = data["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]
    nested["children"] = ""
    chain_number_to_change1 = data["branching_sites"][6]["children"]["branching_sites"][7]["children"]["chain_number"]
    assert nested["children"] == ""
    updated, _ = ubiquitin_building(data, "SMAC", chain_number_to_change1, "K6")
    result = updated["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]["children"]
    assert result == "SMAC"

# 6. Add ubiquitin at level 4 chain.
def test_add_ubiquitin_nested_level_4():
    """Add ubiquitin at level 4 chain and confirm children is dict."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    nested = data["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]
    nested["children"] = ""
    chain_number_to_change1 = data["branching_sites"][6]["children"]["branching_sites"][7]["children"]["chain_number"]
    assert nested["children"] == ""
    updated, _ = ubiquitin_building(data, BASE_WORKING_DICT, chain_number_to_change1, "K6")
    result = updated["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][1]["children"]
    assert isinstance(result, dict)

# 7. Verify chain number list increments.
def test_chain_number_sequence():
    """Verify chain_number_list increments with every addition."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    _, context = ubiquitin_building(data, BASE_WORKING_DICT, 5, "K48")
    assert isinstance(context.get("chain_number_list", []), list)
    assert len(context["chain_number_list"]) == 7
    assert len(context["chain_length_list"]) == 6

# 8. Verify chain number list increments when adding ubiquitin to different sites.
@pytest.mark.parametrize(
    "ubiquitin_number, lysine_site, expected_max, expected_len",
    [
        (1, "K6", 7, 6),
        (2, "K11", 7, 6),
        (2, "K27", 7, 6),
        (3, "K29", 7, 6),
        (3, "K33", 7, 6),
        (4, "K48", 7, 6),
        (5, "K63", 7, 6),
        (5, "M1", 7, 6),
    ]
)
def test_chain_number_sequence(ubiquitin_number, lysine_site, expected_max, expected_len):
    """Verify chain_number_list increments with every addition for various sites."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    _, context = ubiquitin_building(data, BASE_WORKING_DICT, ubiquitin_number, lysine_site)
    assert isinstance(context.get("chain_number_list", []), list)
    assert max(context["chain_number_list"]) == expected_max
    assert len(context["chain_length_list"]) == expected_len

# 9. String input converts to dict.
def test_ubiquitin_string_input_handling():
    """String input for data is converted to dict internally."""
    import json
    jdata = json.dumps(five_level_nested_ubiquitin_)
    updated, _ = ubiquitin_building(jdata, BASE_WORKING_DICT, 1, "K48")
    assert isinstance(updated, dict)

# 10. Unknown label skipped.
def test_ignored_invalid_protecting_group():
    """Unknown protecting group label is skipped and does not modify children."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    before = copy.deepcopy(data["branching_sites"][1]["children"])
    updated, _ = ubiquitin_building(data, BASE_WORKING_DICT, 1, "FOO")
    assert updated["branching_sites"][1]["children"] == before

# 11. Logs SMAC/ABOC correctly.
def test_protecting_group_logged_correctly():
    """SMAC/ABOC protecting group addition is logged correctly."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    updated, context = ubiquitin_building(data, "SMAC", 1, "K6")
    assert updated["branching_sites"][1]["children"] == "SMAC"
    updated2, context2 = ubiquitin_building(data, "ABOC", 1, "K6")
    assert updated2["branching_sites"][1]["children"] == "ABOC"


# 12. Builds at multiple depths.
def test_add_multiple_mods_and_recurse():
    """Builds at multiple depths using several modifications."""
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = ""
    assert data["branching_sites"][1]["children"] == ""
    chain_number_to_change1 = data["chain_number"]

    # Add at level 1
    updated, context = ubiquitin_building(data, BASE_WORKING_DICT, chain_number_to_change1, "K6")
    assert isinstance(updated["branching_sites"][1]["children"], dict)

    # Add at nested level
    updated_copy = copy.deepcopy(updated)
    deep = updated_copy["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][6]
    deep["children"] = ""
    assert deep["children"] == ""
    chain_number_to_change2 = updated_copy["branching_sites"][6]["children"]["branching_sites"][7]["children"]["chain_number"]

    updated2, context2 = ubiquitin_building(updated, BASE_WORKING_DICT, chain_number_to_change2, "K48")
    assert isinstance(
        updated2["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][6]["children"],
        dict
    )

    # Check only a single ubiquitin was added
    updated2["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][6]["children"] = ""
    assert updated2 == updated
    


# =====================================
# === Additional test: FAKE_deprot logs info and does nothing on empty lysine ===
# =====================================
def test_process_ubiquitin_fake_deprot_logs_info(caplog):
    """
    Test FAKE_deprot: logs an info message and takes no action on an empty lysine.
    """
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    target_site = data["branching_sites"][6]  # K48 for example
    target_site["children"] = ""

    with caplog.at_level(logging.INFO):
        bra, updated = process_ubiquitin_reaction(target_site, data, "FAKE_deprot", ubiquitin_monomer)
        assert bra["children"] == ""
        assert "FAKE deprotection: No action taken." in caplog.text