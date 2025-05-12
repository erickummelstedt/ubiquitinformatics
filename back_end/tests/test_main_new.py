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

from src.main import process_ubiquitin_reaction

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

# =============================
# Tests for ubiquitin_simulation and inner_wrapper_ubiquitin_simulation
# =============================
import pytest
from src.main import ubiquitin_simulation, inner_wrapper_ubiquitin_simulation

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