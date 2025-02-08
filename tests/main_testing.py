import pytest
import json
import copy

from copy import deepcopy

import sys

# home_dir = os.path.expanduser('~')
local_path = '/Users/ekummelstedt/le_code_base/ubiquitinformatics'
sys.path.insert(0, local_path)

# Import the functions from the original code
from src.main_testing import find_branching_site, relabelling_ubiquitin_numbers, inner_wrapper_relabelling_ubiquitin_numbers, validate_branching_sites
from src.main_testing import getting_multimer_string_name

# Sample deeply nested ubiquitin dictionary for testing
k48_dimer_ubiquitin = {
    "protein": "1ubq-histag",
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
        {"site_name": "K48","sequence_id": "FAG(K)QLE","children": {"protein": "1ubq",
                                                                    "chain_number": 2,
                                                                    "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                                                    "chain_length": 76,
                                                                    "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                                                                        {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                                                                        {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                                                                        {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                                                                        {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                                                                        {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                                                                        {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                                                                                        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}]}},
        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}]}

string_k48_dimer_ubiquitin = str(k48_dimer_ubiquitin)

ubiquitin_monomer = {
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
                        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}]}

histag_ubiquitin_monomer = {
    "protein": "1ubq",
    "chain_number": 1,
    "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
    "chain_length": 83,
    "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                        {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                        {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                        {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                        {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                        {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                        {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}]}

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

def test_relabelling_ubiquitin_numbers():
    """
    Test `relabelling_ubiquitin_numbers` to ensure proper renumbering and chain 
    length adjustments in a nested ubiquitin structure.
    
    Steps:
    1. Create a deep copy of the `k48_dimer_ubiquitin` dictionary to avoid modifying original data.
    2. Run `relabelling_ubiquitin_numbers` on the test copy.
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
    updated_ubiquitin = relabelling_ubiquitin_numbers(test_ubiquitin)

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
    Test the `relabelling_ubiquitin_numbers` function to ensure correct renumbering of ubiquitin chains
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
    updated_ubiquitin = relabelling_ubiquitin_numbers(test_ubiquitin)
    # Step 4:
    assert updated_ubiquitin["branching_sites"][6]["children"]["branching_sites"][7]["children"]["chain_number"] == 3, \
        "Expected chain number 3 for the first nested ubiquitin monomer"
    assert updated_ubiquitin["branching_sites"][6]["children"]["branching_sites"][7]["children"]["branching_sites"][6]["children"]["chain_number"] == 4, \
        "Expected chain number 4 for the deeply nested ubiquitin monomer"

def test_empty_input():
    """
    Test that `relabelling_ubiquitin_numbers` raises KeyError when given an empty dictionary as input.
    """
    with pytest.raises(KeyError, match="Missing required keys"):
        relabelling_ubiquitin_numbers({})


def test_missing_keys():
    """
    Test that `relabelling_ubiquitin_numbers` raises a KeyError when required keys are missing from the input dictionary.
    """
    incomplete_dict = {"protein": "1ubq"}  # Missing required keys

    with pytest.raises(KeyError, match="Missing required keys"):
        relabelling_ubiquitin_numbers(incomplete_dict)


# Sample ubiquitin data in JSON string format
sample_ubiquitin_json = json.dumps(k48_dimer_ubiquitin)

@pytest.mark.parametrize("input_data, expected_chain_number, expected_length", [
    (k48_dimer_ubiquitin, 1, 83),
    (sample_ubiquitin_json, 1, 83),
])

def test_json_loading_valid(input_data, expected_chain_number, expected_length):
    """
    Test relabelling_ubiquitin_numbers with both dictionary and JSON string input.
    """
    result = relabelling_ubiquitin_numbers(input_data)
    assert result["chain_number"] == expected_chain_number
    assert result["chain_length"] == expected_length


def test_loading_is_string():
    """
    Test relabelling_ubiquitin_numbers with invalid JSON string input.
    Test that the function pulls out the JSON.
    """
    invalid_json = copy.deepcopy(string_k48_dimer_ubiquitin)  # Improper JSON format (single quotes)
    
    result = relabelling_ubiquitin_numbers(invalid_json)
    assert result["chain_number"] == 1
    assert result["chain_length"] == 83

def test_json_loading_empty():
    """
    Test relabelling_ubiquitin_numbers with empty JSON string.
    """
    empty_json = "{}"
    with pytest.raises(KeyError):
        relabelling_ubiquitin_numbers(empty_json)


def test_json_loading_partial_data():
    """
    Test relabelling_ubiquitin_numbers with missing keys in JSON.
    """
    incomplete_json = json.dumps({"protein": "1ubq"})
    with pytest.raises(KeyError):
        relabelling_ubiquitin_numbers(incomplete_json)


def test_json_loading_non_json_string():
    """
    Test relabelling_ubiquitin_numbers with a non-JSON string that is not representative of the correct Dictionary structure.
    """
    non_json_input = "This is not JSON"
    with pytest.raises(ValueError):
        relabelling_ubiquitin_numbers(non_json_input)


def test_json_loading_with_extra_keys():
    """
    Test relabelling_ubiquitin_numbers with additional unexpected keys in JSON. 
    Assert the the additional unexpected keys create errors.
    """
    erroneous_json = json.dumps({
        "protein": "1ubq",
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH",
        "unexpected_key": "some_value"
    })
    
    with pytest.raises(KeyError):
        relabelling_ubiquitin_numbers(erroneous_json)

def test_json_loading_with_unexpected_keys():
    """
    Test relabelling_ubiquitin_numbers with JSON containing unexpected values as keys.
    """
    invalid_json = copy.deepcopy(k48_dimer_ubiquitin)
    invalid_json[123]= "value"
    with pytest.raises(KeyError, match="Unexpected keys found"):
        relabelling_ubiquitin_numbers(invalid_json)

def test_deeply_nested_erroneous_key():
    """
    Test `relabelling_ubiquitin_numbers` to ensure it raises KeyError when an 
    erroneous key exists deep inside the nested ubiquitin structure.
    """
    # Create a deep copy of valid ubiquitin data
    test_ubiquitin = copy.deepcopy(k48_dimer_ubiquitin)

    # Introduce an erroneous key deep in the nested structure
    test_ubiquitin["branching_sites"][6]["children"]['chain_number1'] = test_ubiquitin["branching_sites"][6]["children"]['chain_number']

    with pytest.raises(KeyError, match="Unexpected keys found"):
        relabelling_ubiquitin_numbers(test_ubiquitin)


# Define valid site names
VALID_SITE_NAMES = {"K6", "K11", "K27", "K29", "K33", "K48", "K63", "M1"}

@pytest.mark.parametrize("site_name, sequence_id, children, expected_error", [
    ("K48", "FAG(K)QLE", "", None),  # ✅ Valid case
    ("K63", "NIQ(K)EST", "SMAC", None),  # ✅ Valid case with SMAC
    ("K11", "LTG(K)TIT", {}, None),  # ✅ Valid case with a nested dictionary
    ("K29", "VKA(K)IQD", "ABOC", None),  # ✅ Valid case with ABOC
    ("K27", "ENV(K)AKI", "INVALID_CHILD", "Invalid children format"),  # ❌ Invalid children
    ("K33", "IQD(K)EGI", None, "Invalid children format"),  # ❌ None is not allowed
    ("K99", "XYZ(K)ABC", "", "Invalid site_name"),  # ❌ Invalid site name
    ("K48", "XYZ_K_ABC", "", "Invalid sequence_id format"),  # ❌ Missing parentheses in sequence_id
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
    }, "his-ubi1[]"),

    # Test Case 2: Ubiquitin Dimer with K48 Linkage
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
    }, "his-ubi1[K48_ubi2[]]"),

    # Test Case 3: Ubiquitin Trimer with K48 and K63 Linkages
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
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
            }}
        ]
    }, "his-ubi1[K48_ubi2[K63_ubi3[]]]"),

    # Test Case 4: Ubiquitin with Protecting Groups (SMAC and ABOC)
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
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "SMAC"},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": "ABOC"}
        ]
    }, "his-ubi1[K48_SMAC-K63_ABOC]"),

    # Test Case 5: Deeply Nested Ubiquitin with All Branching Sites Used
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
    }, "his-ubi1[K48_ubi2[K48_ubi3[K48_ubi4[]]-K63_SMAC]]"),
])
def test_getting_multimer_string_name(ubiquitin_structure, expected_multimer_string):
    """
    Test `getting_multimer_string_name` to ensure correct string formation 
    for ubiquitin structures with all possible branching sites.
    """
    result = getting_multimer_string_name(copy.deepcopy(ubiquitin_structure))
    assert result == expected_multimer_string, f"Expected: {expected_multimer_string}, Got: {result}"


## IMPORT ##
from src.main_testing import find_number_of_ABOC_SMAC, find_number_of_ABOC, find_number_of_SMAC

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
def test_find_number_of_ABOC_SMAC(ubiquitin_structure, expected_aboc, expected_smac):
    """
    Test `find_number_of_ABOC_SMAC()` to ensure it correctly counts 
    the number of ABOC and SMAC protecting groups.
    """
    aboc_count, smac_count = find_number_of_ABOC_SMAC(copy.deepcopy(ubiquitin_structure))
    assert len(aboc_count) == expected_aboc, f"Expected {expected_aboc} ABOC, got {len(aboc_count)}"
    assert len(smac_count) == expected_smac, f"Expected {expected_smac} SMAC, got {len(smac_count)}"




import pytest
import copy
from src.main_testing import find_max_chain_number

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
def test_find_max_chain_number(ubiquitin_structure, expected_max_chain):
    assert find_max_chain_number(copy.deepcopy(ubiquitin_structure)) == expected_max_chain

from src.main_testing import validate_all_branching_sites

@pytest.mark.parametrize("ubiquitin_structure, should_raise, expected_missing_sites, expected_chain_number", [
    # ✅ Test Case 1: Valid Monomer (No missing sites)
    ({
        "protein": "1ubq",
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
    }, False, None, None),

    # ❌ Test Case 2: Missing K29 in a Monomer
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, True, {"K29"}, 1),

    # ❌ Test Case 3: Missing M1 & K29 in a Nested Dimer
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "branching_sites": [
                    {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
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
    }, True, {"M1", "K29"}, 2),

    # ❌ Test Case 4: Deeply Nested Ubiquitin Missing K48
    ({
        "protein": "1ubq",
        "chain_number": 1,
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
                        "branching_sites": [
                            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                        ]
                    }},
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, True, {"K48"}, 3),

    # ✅ Test Case 5: Correctly Nested Trimer (All Sites Present)
    ({
        "protein": "1ubq",
        "chain_number": 1,
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
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}

        ]
    }, False, None, None),
])
def test_validate_all_branching_sites(ubiquitin_structure, should_raise, expected_missing_sites, expected_chain_number):
    """
    Tests validate_all_branching_sites() to ensure that missing sites are detected at all levels.
    """
    if should_raise:
        with pytest.raises(AssertionError, match=f"Missing sites {expected_missing_sites} in Ubiquitin {expected_chain_number}"):
            validate_all_branching_sites(copy.deepcopy(ubiquitin_structure))
    else:
        validate_all_branching_sites(copy.deepcopy(ubiquitin_structure))  # Should not raise any error