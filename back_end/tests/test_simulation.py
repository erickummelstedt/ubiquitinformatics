import types
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
    find_max_chain_number, \
    process_ubiquitin_reaction, \
    ubiquitin_simulation, \
    inner_wrapper_ubiquitin_simulation, \
    handle_lysine_modification, \
    ubiquitin_building, \
    inner_wrapper_ubiquitin_building

from src.simulation import simulate_E2_step, simulate_deprot_step

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

# ------------------- Simulate Reactions Step Tests -------------------
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


# ------------------- Simulate Deprot Step Tests -------------------

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
