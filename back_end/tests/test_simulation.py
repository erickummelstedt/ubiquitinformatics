import types
import pytest
import json
import copy
from copy import deepcopy
import sys
import pytest

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
    process_ubiquitin_reaction, \
    ubiquitin_simulation, \
    inner_wrapper_ubiquitin_simulation, \
    handle_lysine_modification, \
    ubiquitin_building, \
    inner_wrapper_ubiquitin_building

from src.simulation import \
    simulate_E2_step, \
    simulate_deprot_step, \
    assign_enzyme, \
    determine_elongation_or_branching, \
    determine_reaction_type, \
    validate_conjugated_lysines, \
    assign_correct_E2_enzyme

from src.utils import \
    match_assertion_error_contains,\
    all_strings_exist, \
    all_strings_exist_in_list, \
    inject_fasta_sequence_at_chain,\
    inject_protein_key,\
    inject_branching_sites,\
    convert_json_to_dict

from tests.test_data import \
    five_level_nested_ubiquitin_,\
    k48_dimer_ubiquitin,\
    string_k48_dimer_ubiquitin,\
    ubiquitin_monomer, \
    histag_ubiquitin_monomer,\
    BASE_WORKING_DICT, \
    BASE_CONTEXT, \
    ubi_ubq_1_K48_SMAC,\
    ubi_ubq_1_K63_SMAC,\
    ubi_ubq_1_K48_SMAC_K63_ABOC,\
    ubi_ubq_1_K48_ABOC_K63_SMAC,\
    ubi_ubq_1_K48_ABOC_K63_ABOC,\
    histag_ubi_ubq_1,\
    histag_ubi_ubq_1_K48_aboc,\
    histag_ubi_ubq_1_K63_aboc

# =====================================
# ======= simulate_E2_step ========
# =====================================
def test_simulate_E2_step_outputs_structure():
    """Test that simulate_E2_step returns correctly structured output."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [
        ubi_ubq_1_K48_SMAC,
        ubi_ubq_1_K63_SMAC
    ]

    result = simulate_E2_step(history_dict, donor_list)

    assert len(result) == 4
    assert all(len(r['ubiquitin_history']) == 2 for r in result)
    assert all(len(r['reaction_history']) == 2 for r in result)
    assert all(len(r['donor_history']) == 2 for r in result)
    assert all(len(r['context_history']) == 2 for r in result)

def test_simulate_E2_step_outputs_structure_5_d():
    """Test that simulate_E2_step returns correctly structured output."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1_K48_aboc],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [
        ubi_ubq_1_K48_SMAC,
        ubi_ubq_1_K63_SMAC,
        ubi_ubq_1_K48_SMAC_K63_ABOC,
        ubi_ubq_1_K48_ABOC_K63_SMAC,
        ubi_ubq_1_K48_ABOC_K63_ABOC
    ]

    result = simulate_E2_step(history_dict, donor_list)

    assert len(result) == 10
    assert all(len(r['ubiquitin_history']) == 2 for r in result)
    assert all(len(r['reaction_history']) == 2 for r in result)
    assert all(len(r['donor_history']) == 2 for r in result)
    assert all(len(r['context_history']) == 2 for r in result)


def test_simulate_E2_step_with_no_monomers():
    """Test behavior when monomer list is empty."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = []

    result = simulate_E2_step(history_dict, donor_list)

    # All output lists should be empty
    assert len(result) == 0


def test_simulate_E2_step_appends_correct_reactions():
    """Test that correct reactions are added to history."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [ubi_ubq_1_K48_SMAC] 

    result = simulate_E2_step(history_dict, donor_list)

    expected = [['', 'K48'], 
                ['', 'K63']]
    assert [r['reaction_history'] for r in result] == expected


def test_simulate_E2_step_multimer_structure():
    """Ensure new_multimer structure includes expected fields."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [ubi_ubq_1_K48_SMAC]
    result = simulate_E2_step(history_dict, donor_list)
    assert 'protein' in result[0]['ubiquitin_history'][-1]
    assert 'branching_sites' in result[0]['ubiquitin_history'][-1]


def test_simulate_E2_step_context_contains_chain_data():
    """Ensure context contains chain_number_list and chain_length_list."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [ubi_ubq_1_K48_SMAC]
    result = simulate_E2_step(history_dict, donor_list)
    context = result[0]['context_history'][-1]
    assert 'chain_number_list' in context
    assert 'chain_length_list' in context


def test_simulate_E2_step_missing_keys_raises_keyerror():
    """Expect KeyError if history_dict missing required keys."""
    history_dict = {
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [ubi_ubq_1_K48_SMAC]
    with pytest.raises(KeyError):
        simulate_E2_step(history_dict, donor_list)


def test_simulate_E2_step_with_unreactive_donor():
    """Ensure function still adds unreactive donors to history."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1_K48_aboc],  # already protected, so no K48 reaction expected
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [ubi_ubq_1_K48_SMAC]
    result = simulate_E2_step(history_dict, donor_list)
    assert len(result) == 2  # Still attempts both K48 and K63
    for r in result:
        assert len(r['ubiquitin_history']) == 2


def test_simulate_E2_step_with_invalid_donor_type():
    """Expect convert_json_to_dict to raise ValueError on invalid donor format."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = ['not a valid JSON']  # Improper donor format
    with pytest.raises(Exception):
        simulate_E2_step(history_dict, donor_list)


def test_simulate_E2_step_duplicate_multimers_allowed():
    """Ensure duplicate multimers are added; no deduplication."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [ubi_ubq_1_K48_SMAC, ubi_ubq_1_K48_SMAC]
    result = simulate_E2_step(history_dict, donor_list)
    assert len(result) == 4  # 2 donors * 2 reactions

# =====================================
# ======= simulate_deprot_step ========
# =====================================

def test_simulate_deprot_step_output_structure():
    """Test that simulate_deprot_step returns correct structure."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1_K48_aboc],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    result = simulate_deprot_step(history_dict)

    assert len(result) == 2
    assert all(len(r['ubiquitin_history']) == 2 for r in result)
    assert all(len(r['reaction_history']) == 2 for r in result)
    assert all(len(r['donor_history']) == 2 for r in result)
    assert all(len(r['context_history']) == 2 for r in result)

def test_simulate_deprot_step_appends_correct_reactions():
    """Check that SMAC_deprot and FAKE_deprot are appended correctly."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1_K48_aboc],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    result = simulate_deprot_step(history_dict)
    assert [r['reaction_history'][-1] for r in result] == ['SMAC_deprot', 'FAKE_deprot']

def test_simulate_deprot_step_multimer_contains_expected_keys():
    """Verify new multimers contain expected structure."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1_K63_aboc],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    result = simulate_deprot_step(history_dict)
    for r in result:
        assert 'protein' in r['ubiquitin_history'][-1]
        assert 'branching_sites' in r['ubiquitin_history'][-1]

def test_simulate_deprot_step_context_structure():
    """Ensure context entries include necessary fields."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    result = simulate_deprot_step(history_dict)
    for r in result:
        assert 'chain_number_list' in r['context_history'][-1]
        assert 'chain_length_list' in r['context_history'][-1]

def test_simulate_deprot_step_missing_keys_raises_error():
    """Check missing keys raise KeyError."""
    history_dict = {
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    with pytest.raises(KeyError):
        simulate_deprot_step(history_dict)


# ------------------- Combined E2 + Deprot Step Test -------------------
def test_simulate_e2_then_deprot_combination():
    """Test combined simulation: E2 step followed by deprot step for each result."""
    history_dict = {
        'ubiquitin_history': [histag_ubi_ubq_1],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [ubi_ubq_1_K48_SMAC, ubi_ubq_1_K63_SMAC]
    e2_results = simulate_E2_step(history_dict, donor_list)

    all_deprot_results = []
    for intermediate in e2_results:
        deprot_results = simulate_deprot_step(intermediate)
        all_deprot_results.extend(deprot_results)

    # Expect 4 E2 results x 2 deprot each = 8 total
    assert len(all_deprot_results) == 8
    for r in all_deprot_results:
        assert len(r['ubiquitin_history']) == 3
        assert len(r['reaction_history']) == 3
        assert len(r['donor_history']) == 3
        assert len(r['context_history']) == 3

# =====================================
# Tests for assign_enzyme
# =====================================

def test_assign_enzyme_k48_elongation():
    """Test that K48 elongation assigns gp78/Ube2g2 enzyme."""
    result = assign_enzyme("K48_reaction", "elongation")
    assert result == "gp78/Ube2g2"


def test_assign_enzyme_k48_branching():
    """Test that K48 branching assigns Ube2K enzyme."""
    result = assign_enzyme("K48_reaction", "branching")
    assert result == "Ube2K"


def test_assign_enzyme_k63_elongation():
    """Test that K63 elongation assigns Ubc13/Mms2 enzyme."""
    result = assign_enzyme("K63_reaction", "elongation")
    assert result == "Ubc13/Mms2"


def test_assign_enzyme_k63_branching():
    """Test that K63 branching also assigns Ubc13/Mms2 enzyme."""
    result = assign_enzyme("K63_reaction", "branching")
    assert result == "Ubc13/Mms2"


def test_assign_enzyme_invalid_combination_raises():
    """Test that an invalid combination raises a TypeError."""
    with pytest.raises(TypeError) as exc_info:
        assign_enzyme("K48_reaction", "invalid_linkage")
    assert "Invalid combination of reaction" in str(exc_info.value)

# ===========================================
# Tests for determine_elongation_or_branching
# ===========================================

# Additional tests for determine_elongation_or_branching
def test_determine_elongation_or_branching_elongation():
    """Test when new bound lysine appears once: should return 'elongation'."""
    product_conjugated_lysines = [
        [1, "K48"], [2, "K63"], [3, "K48"]
    ]
    new_bound_lysine = [3, "K48"]
    product_conjugated_lysines.append([new_bound_lysine])
    result = determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert result == "elongation"


def test_determine_elongation_or_branching_branching():
    """Test when new bound lysine appears twice: should return 'branching'."""
    product_conjugated_lysines = [
        [1, "K48"], [2, "K63"], [2, "K48"]
    ]
    new_bound_lysine = [2, "K48"]
    product_conjugated_lysines.append([new_bound_lysine])
    result = determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert result == "branching"


def test_determine_elongation_or_branching_invalid_count_zero():
    """Test when new bound lysine does not exist at all: should raise TypeError."""
    product_conjugated_lysines = [
        [1, "K48"], [2, "K63"]
    ]
    new_bound_lysine = [3, "K48"]
    product_conjugated_lysines.append([new_bound_lysine])
    with pytest.raises(TypeError) as exc_info:
        determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert "expected 1 or 2" in str(exc_info.value)


def test_determine_elongation_or_branching_invalid_count_more_than_two():
    """Test when new bound lysine appears more than twice: should raise TypeError."""
    product_conjugated_lysines = [
        [1, "K63"], [2, "K63"], [2, "K48"], [2, "K29"]
    ]
    new_bound_lysine = [2, "K63"]
    product_conjugated_lysines.append([new_bound_lysine])
    with pytest.raises(TypeError) as exc_info:
        determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert "expected 1 or 2" in str(exc_info.value)


def test_determine_elongation_or_branching_with_empty_product_list():
    """Test when product_conjugated_lysines is length 1 and new lysine is added once."""
    product_conjugated_lysines = [[1, "K48"]]
    new_bound_lysine = [1, "K48"]
    product_conjugated_lysines.append([new_bound_lysine])
    result = determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert result == "elongation"

# =====================================
# Tests for determine_reaction_type
# =====================================

def test_determine_reaction_type_k48():
    """Test that K48 lysine returns correct reaction type."""
    new_bound_lysine = [3, "K48"]
    result = determine_reaction_type(new_bound_lysine)
    assert result == "K48_reaction"


def test_determine_reaction_type_k63():
    """Test that K63 lysine returns correct reaction type."""
    new_bound_lysine = [5, "K63"]
    result = determine_reaction_type(new_bound_lysine)
    assert result == "K63_reaction"


def test_determine_reaction_type_invalid_lysine_raises():
    """Test that an invalid lysine site raises TypeError."""
    new_bound_lysine = [2, "K33"]
    with pytest.raises(TypeError) as exc_info:
        determine_reaction_type(new_bound_lysine)
    assert "does not contain K48 or K63" in str(exc_info.value)


def test_determine_reaction_type_non_string_site_raises():
    """Test that a non-string lysine site raises TypeError."""
    new_bound_lysine = [4, 48]
    with pytest.raises(TypeError) as exc_info:
        determine_reaction_type(new_bound_lysine)
    assert "does not contain K48 or K63" in str(exc_info.value)

# =====================================
# Tests for validate_conjugated_lysines
# =====================================

def test_validate_conjugated_lysines_all_valid():
    """Test that validation passes when only K48 and K63 are present."""
    context = {
        "conjugated_lysines": [[1, "K48"], [2, "K63"], [3, "K48"]]
    }
    # Should not raise any exceptions
    validate_conjugated_lysines(context)


def test_validate_conjugated_lysines_with_invalid_entry():
    """Test that validation raises TypeError when unsupported lysines are included."""
    context = {
        "conjugated_lysines": [[1, "K48"], [2, "K33"], [3, "K63"]]
    }
    with pytest.raises(TypeError) as exc_info:
        validate_conjugated_lysines(context)
    assert "unsupported conjugated lysines" in str(exc_info.value)


def test_validate_conjugated_lysines_with_empty_list():
    """Test that validation passes when the list is empty."""
    context = {
        "conjugated_lysines": []
    }
    # Should not raise any exceptions
    validate_conjugated_lysines(context)


def test_validate_conjugated_lysines_with_malformed_entries():
    """Test that malformed entries are ignored, and validation passes if remaining are valid."""
    context = {
        "conjugated_lysines": [[1, "K48"], ["invalid"], [2, "K63"], [3], "random"]
    }
    # Should not raise an error because malformed entries are skipped
    validate_conjugated_lysines(context)


# ============================================================
# Build reaction sequence for testing assign_correct_E2_enzyme
# ============================================================
# === First Reaction Path ===

# Dimer formation (K48 elongation)
acceptor_1_0, context_1_0 = iterate_through_ubiquitin(histag_ubi_ubq_1_K63_aboc)
acceptor = acceptor_1_0.copy()
context = context_1_0.copy()
donor = ubi_ubq_1_K48_SMAC
reaction = "K48"
acceptor_1_1a, context_1_1a = ubiquitin_simulation(acceptor, donor, reaction)

# Dimer deprotection
acceptor, context = acceptor_1_1a.copy(), context_1_1a.copy()
acceptor_1_1b, context_1_1b = ubiquitin_simulation(acceptor, '', "FAKE_deprot")

# Trimer formation (K63 elongation)
acceptor, context = acceptor_1_1b.copy(), context_1_1b.copy()
donor = ubi_ubq_1_K48_ABOC_K63_ABOC
reaction = "K63"
acceptor_1_2a, context_1_2a = ubiquitin_simulation(acceptor, donor, reaction)

# Trimer deprotection
acceptor, context = acceptor_1_2a.copy(), context_1_2a.copy()
acceptor_1_2b, context_1_2b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")

# Tetramer formation (K48 branching)
acceptor, context = acceptor_1_2b.copy(), context_1_2b.copy()
donor = ubi_ubq_1_K48_ABOC_K63_ABOC
reaction = "K48"
acceptor_1_3a, context_1_3a = ubiquitin_simulation(acceptor, donor, reaction)

# Tetramer deprotection
acceptor, context = acceptor_1_3a.copy(), context_1_3a.copy()
acceptor_1_3b, context_1_3b = ubiquitin_simulation(acceptor, '', "FAKE_deprot")


# === Second Reaction Path ===

# Dimer formation (K48 elongation)
acceptor_2_0, context_2_0 = iterate_through_ubiquitin(histag_ubi_ubq_1_K63_aboc)
acceptor, context = acceptor_2_0.copy(), context_2_0.copy()
donor = ubi_ubq_1_K48_ABOC_K63_SMAC
reaction = "K48"
acceptor_2_1a, context_2_1a = ubiquitin_simulation(acceptor, donor, reaction)

# Dimer deprotection
acceptor, context = acceptor_2_1a.copy(), context_2_1a.copy()
acceptor_2_1b, context_2_1b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")

# Trimer formation (K63 elongation)
acceptor, context = acceptor_2_1b.copy(), context_2_1b.copy()
donor = ubi_ubq_1_K63_SMAC
reaction = "K63"
acceptor_2_2a, context_2_2a = ubiquitin_simulation(acceptor, donor, reaction)

# Trimer deprotection
acceptor, context = acceptor_2_2a.copy(), context_2_2a.copy()
acceptor_2_2b, context_2_2b = ubiquitin_simulation(acceptor, '', "FAKE_deprot")

# Tetramer formation (K48 branching)
acceptor, context = acceptor_2_2b.copy(), context_2_2b.copy()
donor = ubi_ubq_1_K48_ABOC_K63_ABOC
reaction = "K48"
acceptor_2_3a, context_2_3a = ubiquitin_simulation(acceptor, donor, reaction)

# Tetramer deprotection
acceptor, context = acceptor_2_3a.copy(), context_2_3a.copy()
acceptor_2_3b, context_2_3b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")

# Pentamer formation (K63 branching)
acceptor, context = acceptor_2_3b.copy(), context_2_3b.copy()
donor = ubi_ubq_1_K48_ABOC_K63_ABOC
reaction = "K63"
acceptor_2_4a, context_2_4a = ubiquitin_simulation(acceptor, donor, reaction)

# Pentamer deprotection
acceptor, context = acceptor_2_4a.copy(), context_2_4a.copy()
acceptor_2_4b, context_2_4b = ubiquitin_simulation(acceptor, '', "FAKE_deprot")

# Parameterized test for assign_correct_E2_enzyme
@pytest.mark.parametrize(
    "reactant, product, expected_enzyme, note",
    [
        (context_1_0, context_1_1a, "gp78/Ube2g2", "Dimer formation via K48 results in elongation"),
        (context_1_1b, context_1_2a, "Ubc13/Mms2", "Trimer formation via K63 results in elongation"),
        (context_1_2b, context_1_3a, "Ube2K", "Tetramer formation via K48 results in branching"),
        (context_2_0, context_2_1a, "gp78/Ube2g2", "Dimer formation via K48 results in elongation (2nd path)"),
        (context_2_1b, context_2_2a, "Ubc13/Mms2", "Trimer formation via K63 results in elongation (2nd path)"),
        (context_2_2b, context_2_3a, "gp78/Ube2g2", "Tetramer formation via K48 results in branching (2nd path)"),
        (context_2_3b, context_2_4a, "Ubc13/Mms2", "Pentamer formation via K63 results in branching"),
    ]
)
def test_assign_correct_E2_enzyme_matches_expected(reactant, product, expected_enzyme, note):
    """Parametrized test: Check assign_correct_E2_enzyme against expected enzyme outcomes."""
    result = assign_correct_E2_enzyme(reactant, product)
    assert result == expected_enzyme, f"{note} failed."

# ============================================================
# Normal testing assign_correct_E2_enzyme
# ============================================================

# Edge case: No difference in conjugated_lysines between reactant and product
def test_assign_correct_E2_enzyme_no_new_conjugation_raises():
    """Test that function raises if no new conjugation site is detected."""
    context = {
        "max_chain_number": 1,
        "conjugated_lysines": [[1, "K48"]],
        "chain_number_list": [1],
        "chain_length_list": []
    }
    with pytest.raises(TypeError) as exc_info:
        assign_correct_E2_enzyme(context, context)
    assert "product_max_chain_number: 1 != reactant_max_chain_number + 1: 2" in str(exc_info.value)