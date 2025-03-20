import pytest
import json
import copy
from copy import deepcopy
import sys

# home_dir = os.path.expanduser('~')
local_path = '/Users/ekummelstedt/le_code_base/ubiquitinformatics'
sys.path.insert(0, local_path)

# Import the functions from the original code
# from src.main_testing import relabelling_ubiquitin_numbers, inner_wrapper_relabelling_ubiquitin_numbers
from src.main import \
    iterate_through_ubiquitin, \
    inner_wrapper_iterate_through_ubiquitin, \
    validate_branching_sites, \
    find_branching_site, \
    validate_protein_keys, \
    affirm_branching_sites

'''
Done
- find_branching_site
- validate_protein_keys
- affirm_branching_sites

Build tests for the following
Everything works; now build tests and clean up code for each of the following functions

- validate_branching_sites
- convert_json_to_dict
- process_current_protein
- process_branch
- log_branching_details
- log_end_of_branching
- log_protein_details
- log_end_of_protein
- find_max_chain_number
- iterate_through_ubiquitin (all tests on deeply nested ubiquitins)


'''

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
        {'chain_number', 'FASTA_sequence', 'chain_length'}
    )
])
def test_validate_protein_keys_missing_keys(invalid_dict, missing_keys):
    """
    Test that `validate_protein_keys` raises a KeyError when required keys are missing.
    """
    allowed_keys = {"protein", "chain_number", "FASTA_sequence", "chain_length", "branching_sites"}

    with pytest.raises(KeyError, match=f"Missing required keys: {missing_keys}. Allowed keys: {allowed_keys}"):
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
        {"extra_key"}
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
        {"extra_key_1", "extra_key_2"}
    )
])

def test_validate_protein_keys_invalid_keys(invalid_dict, invalid_keys):
    """
    Test that `validate_protein_keys` raises a KeyError when invalid keys are present.
    """
    allowed_keys = {"protein", "chain_number", "FASTA_sequence", "chain_length", "branching_sites"}
    
    with pytest.raises(KeyError, match=f"Invalid keys found: {invalid_keys}. Allowed keys: {allowed_keys}"):
        validate_protein_keys(invalid_dict)


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


@pytest.mark.parametrize("ubiquitin_structure, should_raise, expected_exception_message", [
    
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
        
    }, False, None),

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
    }, True, "Missing required sites: {'K29', 'K33', 'K48', 'K63'} in Ubiquitin 2"),
    
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
    }, True, "Invalid sites found: {'K99'} in Ubiquitin 3"),

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
            {"site_name": "K100","sequence_id": "NIQ(K)EST","children": ""}]  # Invalid site
    }, True, "Missing required sites: {'K27', 'K29', 'K48', 'K63'}. Invalid sites found: {'K100'} in Ubiquitin 4"),

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
    }, False, None)

])
def test_affirm_branching_sites(ubiquitin_structure, should_raise, expected_exception_message):
    """
    Test `affirm_branching_sites()` with various ubiquitin structures to verify correct
    validation of required branching sites.
    """
    if should_raise:
        with pytest.raises(KeyError, match=expected_exception_message):
            affirm_branching_sites(ubiquitin_structure)
    else:
        # Should not raise an error
        affirm_branching_sites(ubiquitin_structure)











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
    }, "GG-(ubi1)"),
    
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
    }, "his-(ubi1)"),

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
    }, "GG-(ubi1<K48_(ubi2)>)"),

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
    }, "his-(ubi1<K48_(ubi2<K63_(ubi3)>)>)"),

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
    }, "his-(ubi1<K48_SMAC><K63_ABOC>)"),

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
    }, "his-(ubi1<K48_(ubi2<K48_(ubi3<K48_(ubi4)>)><K63_SMAC>)>)"),
])
def test_iterate_through_ubiquitin(ubiquitin_structure, expected_multimer_string):
    """
    Test `inner_wrapper_iterate_through_ubiquitin` to ensure correct string formation 
    for ubiquitin structures with all possible branching sites.
    """
    result, context = iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))
    multimer_string = context["multimer_string_name"]
    assert multimer_string == expected_multimer_string, f"Expected: {expected_multimer_string}, Got: {result}"





## IMPORT ##
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


import pytest
import copy
from src.main import find_max_chain_number

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

    updated_ubiquitin, context = iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))
    
    assert find_max_chain_number(context) == expected_max_chain

from src.main import affirm_branching_sites

@pytest.mark.parametrize("ubiquitin_structure, should_raise, expected_missing_sites, expected_chain_number", [
    # ✅ Test Case 1: Valid Monomer (No missing sites)
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
    }, False, None, None),

    # ❌ Test Case 2: Missing K29 in a Monomer
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
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
        ]
    }, True, {"K29"}, 1),

    # ❌ Test Case 3: Missing M1 & K29 in a Nested Dimer
    ({
        "protein": "1ubq",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 76,
        "branching_sites": [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": {
                "protein": "1ubq",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                "chain_length": 76,
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
    }, True, {"K29", "M1"}, 2),

    # ❌ Test Case 4: Deeply Nested Ubiquitin Missing K48
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
                    {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}

                ]
            }},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}

        ]
    }, False, None, None),
])
def test_affirm_branching_sites(ubiquitin_structure, should_raise, expected_missing_sites, expected_chain_number):
    """
    Tests affirm_branching_sites() within iterate_through_ubiquitin to ensure that missing sites are detected at all levels.
    """
    allowed_sites = {"M1", "K6", "K11", "K27", "K29", "K33", "K48", "K63"}

    if should_raise:
        with pytest.raises(KeyError, match=f"Missing required sites: {expected_missing_sites} in Ubiquitin {expected_chain_number}. Allowed sites: {allowed_sites}."):
            iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))
    else:
        iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))  # Should not raise any error


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
    all conjugated lysines from ubiquitin structures.
    """
    _, context = iterate_through_ubiquitin(copy.deepcopy(ubiquitin_structure))
    conjugated_lysines = context['conjugated_lysines']

    assert sorted(conjugated_lysines) == sorted(expected_conjugated_lysines), \
        f"Expected {expected_conjugated_lysines}, got {conjugated_lysines}"


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