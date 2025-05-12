import pytest
import json
import copy
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
    find_max_chain_number

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

# could call src.main == src.modification instead

# =============================
# Tests for K_residue_ubi_addition (from src.main)
# =============================

from src.main import K_residue_ubi_addition
from src.utils import convert_json_to_dict
import copy

# Test that a ubiquitin is added to a lysine with no existing children
def test_add_to_unoccupied_site():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    updated = K_residue_ubi_addition(data, 1, "IFV(K)TLT", BASE_WORKING_DICT)
    assert isinstance(updated["branching_sites"][1]["children"], dict)

# Test that a ubiquitin is added even if the site has SMAC or ABOC as a string
def test_add_to_site_with_string_children():
    for protecting_group in ["SMAC", "ABOC"]:
        data = copy.deepcopy(five_level_nested_ubiquitin_)
        data["branching_sites"][1]["children"] = protecting_group
        updated = K_residue_ubi_addition(data, 1, "IFV(K)TLT", BASE_WORKING_DICT)
        assert isinstance(updated["branching_sites"][1]["children"], dict)

# Test that no update occurs if the sequence_id doesn't match any lysine
def test_no_change_if_sequence_id_mismatch():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    original = copy.deepcopy(data["branching_sites"])
    updated = K_residue_ubi_addition(data, 1, "XYZ(K)ABC", BASE_WORKING_DICT)
    assert updated["branching_sites"] == original

# Test that an error is raised if the site already contains a conjugated ubiquitin (dict)
def test_error_on_already_conjugated_site():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    data["branching_sites"][1]["children"] = {"mock": "ubi"}
    with pytest.raises(TypeError, match="already conjugated"):
        K_residue_ubi_addition(data, 1, "IFV(K)TLT", BASE_WORKING_DICT)

# Test that string input in JSON format is correctly converted to a dictionary
def test_json_input_converted():
    import json
    json_data = json.dumps(five_level_nested_ubiquitin_)
    updated = K_residue_ubi_addition(json_data, 1, "IFV(K)TLT", BASE_WORKING_DICT)
    assert isinstance(updated, dict)

# Test that a TypeError is raised if the FASTA sequence does not end in RLRGG
def test_addition_raises_if_RLRGG_tail_missing():
    non_matching_tail = copy.deepcopy(BASE_WORKING_DICT)
    non_matching_tail["FASTA_sequence"] = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGXXXXXX"
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    with pytest.raises(TypeError, match="does not end with RLRGG"):
        K_residue_ubi_addition(data, 1, "IFV(K)TLT", non_matching_tail)

# Test that the correct index is targeted and updated within branching_sites
def test_loop_index_resolved_correctly():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    # Ensure the correct index (1) is updated
    updated = K_residue_ubi_addition(data, 1, "IFV(K)TLT", BASE_WORKING_DICT)
    assert updated["branching_sites"][1]["site_name"] == "K6"
    assert isinstance(updated["branching_sites"][1]["children"], dict)

# Test that a ubiquitin can be added to a deeply nested chain structure
def test_addition_on_deep_nested_chain():
    data = copy.deepcopy(five_level_nested_ubiquitin_)
    deep_branch = data["branching_sites"][6]["children"]["branching_sites"][7]["children"]
    assert deep_branch["chain_number"] == 3
    modified = K_residue_ubi_addition(deep_branch, 3, "ENV(K)AKI", BASE_WORKING_DICT)
    assert isinstance(modified["branching_sites"][3]["children"], dict)

# Test that a TypeError is raised if the input is not a dictionary or JSON string
def test_addition_raises_if_not_dict_or_json():
    invalid_input = "This is not a valid input"
    with pytest.raises(ValueError, match="Invalid JSON format: Unable to parse the string"):
        K_residue_ubi_addition(invalid_input, 1, "IFV(K)TLT", BASE_WORKING_DICT)

# Test that a TypeError is raised if the input is an empty string
def test_addition_raises_if_empty_string():
    empty_input = ""
    with pytest.raises(ValueError, match="Invalid JSON format: Unable to parse the string"):
        K_residue_ubi_addition(empty_input, 1, "IFV(K)TLT", BASE_WORKING_DICT)

# Test that a TypeError is raised if the input is None
def test_addition_raises_if_none():
    none_input = None
    with pytest.raises(TypeError, match="Input must be a dictionary or a JSON string"):
        K_residue_ubi_addition(none_input, 1, "IFV(K)TLT", BASE_WORKING_DICT)

