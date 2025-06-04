import types
import pytest
import json
import copy
from copy import deepcopy
import sys
import pytest
from pathlib import Path

# Dynamically get the backend path relative to this file
current_file = Path(__file__).resolve()
project_root = current_file.parents[2]  # Go up to project root
sys.path.insert(0, str(project_root))
local_path = project_root / 'back_end'

sys.path.insert(0, str(local_path))

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
    simulate_E2_steps, \
    simulate_deprot_steps, \
    assign_enzyme, \
    determine_elongation_or_branching, \
    determine_reaction_type, \
    validate_conjugated_lysines, \
    assign_correct_E2_enzyme, \
    determine_Ube2K_elongation

from src.utils.utils import \
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

# ==========================================================
# === Complex Reaction Path 1 for testing with histories ===
# ==========================================================

acceptor_3_0, context_3_0 = iterate_through_ubiquitin(histag_ubi_ubq_1_K63_aboc)
history_dict_3_0 = {}
history_dict_3_0['ubiquitin_history'] = [acceptor_3_0]
history_dict_3_0['reaction_history'] = ['']
history_dict_3_0['donor_history'] = ['']
history_dict_3_0['context_history'] = [context_3_0]

# Dimer formation (K48 elongation)
acceptor, context = acceptor_3_0.copy(), context_3_0.copy()
donor = ubi_ubq_1_K48_ABOC_K63_SMAC
reaction = "K48"
acceptor_3_1a, context_3_1a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_3_1a = {}
history_dict_3_1a['ubiquitin_history'] = history_dict_3_0['ubiquitin_history'] + [acceptor_3_0]
history_dict_3_1a['reaction_history'] = history_dict_3_0['reaction_history'] + [reaction]
history_dict_3_1a['donor_history'] = history_dict_3_0['donor_history'] + [donor]
history_dict_3_1a['context_history'] = history_dict_3_0['context_history'] + [context_3_0]


# Dimer deprotection
acceptor, context = acceptor_3_1a.copy(), context_3_1a.copy()
acceptor_3_1b, context_3_1b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")
history_dict_3_1b = {}
history_dict_3_1b['ubiquitin_history'] = history_dict_3_1a['ubiquitin_history'] + [acceptor_3_1b]
history_dict_3_1b['reaction_history'] = history_dict_3_1a['reaction_history'] + [reaction]
history_dict_3_1b['donor_history'] = history_dict_3_1a['donor_history'] + [donor]
history_dict_3_1b['context_history'] = history_dict_3_1a['context_history'] + [context_3_1b]


# Trimer formation (K63 elongation)
acceptor, context = acceptor_3_1b.copy(), context_3_1b.copy()
donor = ubi_ubq_1_K63_SMAC
reaction = "K63"
acceptor_3_2a, context_3_2a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_3_2a = {}
history_dict_3_2a['ubiquitin_history'] = history_dict_3_1b['ubiquitin_history'] + [acceptor_3_2a]
history_dict_3_2a['reaction_history'] = history_dict_3_1b['reaction_history'] + [reaction]
history_dict_3_2a['donor_history'] = history_dict_3_1b['donor_history'] + [donor]
history_dict_3_2a['context_history'] = history_dict_3_1b['context_history'] + [context_3_2a]


# Trimer deprotection
acceptor, context = acceptor_3_2a.copy(), context_3_2a.copy()
acceptor_3_2b, context_3_2b = ubiquitin_simulation(acceptor, '', "FAKE_deprot")
history_dict_3_2b = {}
history_dict_3_2b['ubiquitin_history'] = history_dict_3_2a['ubiquitin_history'] + [acceptor_3_2b]
history_dict_3_2b['reaction_history'] = history_dict_3_2a['reaction_history'] + [reaction]
history_dict_3_2b['donor_history'] = history_dict_3_2a['donor_history'] + [donor]
history_dict_3_2b['context_history'] = history_dict_3_2a['context_history'] + [context_3_2b]


# Tetramer formation (K48 branching)
acceptor, context = acceptor_3_2b.copy(), context_3_2b.copy()
donor = ubi_ubq_1_K48_ABOC_K63_ABOC
reaction = "K48"
acceptor_3_3a, context_3_3a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_3_3a = {}
history_dict_3_3a['ubiquitin_history'] = history_dict_3_2b['ubiquitin_history'] + [acceptor_3_3a]
history_dict_3_3a['reaction_history'] = history_dict_3_2b['reaction_history'] + [reaction]
history_dict_3_3a['donor_history'] = history_dict_3_2b['donor_history'] + [donor]
history_dict_3_3a['context_history'] = history_dict_3_2b['context_history'] + [context_3_3a]


# Tetramer deprotection
acceptor, context = acceptor_3_3a.copy(), context_3_3a.copy()
acceptor_3_3b, context_3_3b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")
history_dict_3_3b = {}
history_dict_3_3b['ubiquitin_history'] = history_dict_3_3a['ubiquitin_history'] + [acceptor_3_3b]
history_dict_3_3b['reaction_history'] = history_dict_3_3a['reaction_history'] + [reaction]
history_dict_3_3b['donor_history'] = history_dict_3_3a['donor_history'] + [donor]
history_dict_3_3b['context_history'] = history_dict_3_3a['context_history'] + [context_3_3b]

# =================================================================================
# === Complex Reaction Path 2 for testing with histories ==========================
# === Two ubiquitins can be added to K48's if you simulate after last trimer deprot
# =================================================================================

acceptor_4_0, context_4_0 = iterate_through_ubiquitin(histag_ubi_ubq_1_K63_aboc)
history_dict_4_0 = {}
history_dict_4_0['ubiquitin_history'] = [acceptor_4_0]
history_dict_4_0['reaction_history'] = ['']
history_dict_4_0['donor_history'] = ['']
history_dict_4_0['context_history'] = [context_4_0]

# Dimer formation (K48 elongation)
acceptor, context = acceptor_4_0.copy(), context_4_0.copy()
donor = ubi_ubq_1_K63_SMAC
reaction = "K48"
acceptor_4_1a, context_4_1a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_4_1a = {}
history_dict_4_1a['ubiquitin_history'] = history_dict_4_0['ubiquitin_history'] + [acceptor_4_0]
history_dict_4_1a['reaction_history'] = history_dict_4_0['reaction_history'] + [reaction]
history_dict_4_1a['donor_history'] = history_dict_4_0['donor_history'] + [donor]
history_dict_4_1a['context_history'] = history_dict_4_0['context_history'] + [context_4_0]


# Dimer deprotection
acceptor, context = acceptor_4_1a.copy(), context_4_1a.copy()
acceptor_4_1b, context_4_1b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")
history_dict_4_1b = {}
history_dict_4_1b['ubiquitin_history'] = history_dict_4_1a['ubiquitin_history'] + [acceptor_4_1b]
history_dict_4_1b['reaction_history'] = history_dict_4_1a['reaction_history'] + [reaction]
history_dict_4_1b['donor_history'] = history_dict_4_1a['donor_history'] + [donor]
history_dict_4_1b['context_history'] = history_dict_4_1a['context_history'] + [context_4_1b]


# Trimer formation (K63 elongation)
acceptor, context = acceptor_4_1b.copy(), context_4_1b.copy()
donor = ubi_ubq_1_K63_SMAC
reaction = "K63"
acceptor_4_2a, context_4_2a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_4_2a = {}
history_dict_4_2a['ubiquitin_history'] = history_dict_4_1b['ubiquitin_history'] + [acceptor_4_2a]
history_dict_4_2a['reaction_history'] = history_dict_4_1b['reaction_history'] + [reaction]
history_dict_4_2a['donor_history'] = history_dict_4_1b['donor_history'] + [donor]
history_dict_4_2a['context_history'] = history_dict_4_1b['context_history'] + [context_4_2a]

# Trimer deprotection
acceptor, context = acceptor_4_2a.copy(), context_4_2a.copy()
acceptor_4_2b, context_4_2b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")
history_dict_4_2b = {}
history_dict_4_2b['ubiquitin_history'] = history_dict_4_2a['ubiquitin_history'] + [acceptor_4_2b]
history_dict_4_2b['reaction_history'] = history_dict_4_2a['reaction_history'] + [reaction]
history_dict_4_2b['donor_history'] = history_dict_4_2a['donor_history'] + [donor]
history_dict_4_2b['context_history'] = history_dict_4_2a['context_history'] + [context_4_2b]

# =====================================
# ======= simulate_E2_steps ========
# =====================================
def test_simulate_E2_steps_outputs_structure():
    """Test that simulate_E2_steps returns correctly structured output."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [
        ubi_ubq_1_K48_SMAC,
        ubi_ubq_1_K63_SMAC
    ]

    result = simulate_E2_steps(history_dict, donor_list)

    assert len(result) == 2
    assert all(len(r['ubiquitin_history']) == 2 for r in result)
    assert all(len(r['reaction_history']) == 2 for r in result)
    assert all(len(r['donor_history']) == 2 for r in result)
    assert all(len(r['context_history']) == 2 for r in result)

def test_simulate_E2_steps_outputs_structure_5_d():
    """
    
    Test that simulate_E2_steps returns correctly structured output.
    
    K48 is blocked so there should only be K63 reactions and ubi_ubq_1_K48_SMAC reacts with itself.
    Thus there are only 4 reactions.
    """
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1_K48_aboc)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [
        ubi_ubq_1_K48_SMAC,
        ubi_ubq_1_K63_SMAC,
        ubi_ubq_1_K48_SMAC_K63_ABOC,
        ubi_ubq_1_K48_ABOC_K63_SMAC,
        ubi_ubq_1_K48_ABOC_K63_ABOC
    ]

    result = simulate_E2_steps(history_dict, donor_list)

    assert len(result) == 4
    assert all(len(r['ubiquitin_history']) == 2 for r in result)
    assert all(len(r['reaction_history']) == 2 for r in result)
    assert all(len(r['donor_history']) == 2 for r in result)
    assert all(len(r['context_history']) == 2 for r in result)


# ------------------- Additional tests for simulate_E2_steps SMAC/ABOC logic -------------------
def test_simulate_E2_steps_skips_K63_reaction_for_K48_SMAC():
    """Test that simulate_E2_steps skips K63 reaction for K48_SMAC monomer."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1_K63_aboc)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [ubi_ubq_1_K48_SMAC]

    result = simulate_E2_steps(history_dict, donor_list)

    assert len(result) == 1
    assert all(len(r['ubiquitin_history']) == 2 for r in result)
    assert all(len(r['reaction_history']) == 2 for r in result)
    assert all(len(r['donor_history']) == 2 for r in result)
    assert all(len(r['context_history']) == 2 for r in result)

def test_simulate_E2_steps_skips_K48_reactions_for_nested_acceptor_with_two_free_K48():
    """
    
    Test that simulate_E2_steps skips reactions that adds two ubiquitin onto an acceptor than 
    has two free K48's.
    
    Additionally tests that simulate_E2_steps skips K63 reaction for K48_SMAC monomer
    
    """
    donor_list = [
        ubi_ubq_1_K48_SMAC,
        ubi_ubq_1_K63_SMAC,
        ubi_ubq_1_K48_SMAC_K63_ABOC,
        ubi_ubq_1_K48_ABOC_K63_SMAC,
        ubi_ubq_1_K48_ABOC_K63_ABOC
    ]

    result = simulate_E2_steps(history_dict_4_2b, donor_list)

    assert len(result) == 4
    assert all(len(r['ubiquitin_history']) == 6 for r in result)
    assert all(len(r['reaction_history']) == 6 for r in result)
    assert all(len(r['donor_history']) == 6 for r in result)
    assert all(len(r['context_history']) == 6 for r in result)

def test_simulate_E2_steps_skips_K48_reaction_for_K63_SMAC():
    """Test that simulate_E2_steps skips K48 reaction for K63_SMAC monomer."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1_K48_aboc)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [ubi_ubq_1_K63_SMAC]
    result = simulate_E2_steps(history_dict, donor_list)

    # Only one valid reaction (K63), so only 1 result expected
    assert len(result) == 1
    assert result[0]['reaction_history'][-1] == 'Ubc13/Mms2'

def test_simulate_E2_steps_both_smac_monomers_valid_reactions():
    """Test that both K48_SMAC and K63_SMAC each trigger one valid reaction."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [ubi_ubq_1_K48_SMAC, ubi_ubq_1_K63_SMAC]
    result = simulate_E2_steps(history_dict, donor_list)

    # Expect 2 total reactions: K48 for K48_SMAC, K63 for K63_SMAC
    assert len(result) == 2
    assert sorted([r['reaction_history'][-1] for r in result]) == ['Ubc13/Mms2', 'gp78/Ube2g2']

def test_simulate_E2_steps_mixed_valid_and_invalid_combinations():
    """Test that simulate_E2_steps handles a mix of valid and invalid reaction-monomer pairs."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [
        ubi_ubq_1_K48_SMAC,  # K48 valid only
        ubi_ubq_1_K63_SMAC,  # K63 valid only
        ubi_ubq_1_K48_ABOC_K63_SMAC  # both valid
    ]
    result = simulate_E2_steps(history_dict, donor_list)

    # Expect 4 total: 1 from first, 1 from second, 2 from third
    assert len(result) == 4
    reactions = [r['reaction_history'][-1] for r in result]
    assert 'gp78/Ube2g2' in reactions
    assert 'Ubc13/Mms2' in reactions
    assert 'Ube2K' not in reactions


def test_simulate_E2_steps_with_no_monomers():
    """Test behavior when monomer list is empty."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = []

    result = simulate_E2_steps(history_dict, donor_list)

    # All output lists should be empty
    assert len(result) == 0


def test_simulate_E2_steps_appends_correct_reactions():
    """Test that correct reactions are added to history."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [ubi_ubq_1_K48_SMAC] 

    result = simulate_E2_steps(history_dict, donor_list)

    expected = [['', 'gp78/Ube2g2']]
    assert [r['reaction_history'] for r in result] == expected


def test_simulate_E2_steps_multimer_structure():
    """Ensure new_multimer structure includes expected fields."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [ubi_ubq_1_K48_SMAC]
    result = simulate_E2_steps(history_dict, donor_list)
    assert 'protein' in result[0]['ubiquitin_history'][-1]
    assert 'branching_sites' in result[0]['ubiquitin_history'][-1]


def test_simulate_E2_steps_context_contains_chain_data():
    """Ensure context contains chain_number_list and chain_length_list."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [ubi_ubq_1_K48_SMAC]
    result = simulate_E2_steps(history_dict, donor_list)
    context = result[0]['context_history'][-1]
    assert 'chain_number_list' in context
    assert 'chain_length_list' in context


def test_simulate_E2_steps_missing_keys_raises_keyerror():
    """Expect KeyError if history_dict missing required keys."""
    history_dict = {
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    donor_list = [ubi_ubq_1_K48_SMAC]
    with pytest.raises(KeyError):
        simulate_E2_steps(history_dict, donor_list)


def test_simulate_E2_steps_with_unreactive_donor():
    """Ensure function still adds doesn't add unreactive donors to history."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1) # already protected, so no K48 reaction expected
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [ubi_ubq_1_K48_SMAC]
    result = simulate_E2_steps(history_dict, donor_list)
    assert len(result) == 1  # Still attempts both K48 and K63
    for r in result:
        assert len(r['ubiquitin_history']) == 2


def test_simulate_E2_steps_with_invalid_donor_type():
    """Expect convert_json_to_dict to raise ValueError on invalid donor format."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = ['not a valid JSON']  # Improper donor format
    with pytest.raises(Exception):
        simulate_E2_steps(history_dict, donor_list)


def test_simulate_E2_steps_duplicate_multimers_allowed():
    """Ensure duplicate multimers are added; no deduplication."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [ubi_ubq_1_K48_SMAC, ubi_ubq_1_K48_SMAC]
    result = simulate_E2_steps(history_dict, donor_list)
    assert len(result) == 2  # 2 donors * 1 reactions (one reaction allowed)

# =====================================
# ======= simulate_deprot_steps ========
# =====================================

def test_simulate_deprot_steps_output_structure():
    """
    Test that simulate_deprot_steps returns correct structure. 
    
    If there is there is not SMAC to deprotect and the reaction is a SMAC_deprot,
    this history_dict should be skipped. Hence there is only one reaction appended here (FAKE_deprot)

    """
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1_K48_aboc)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    result = simulate_deprot_steps(history_dict)

    assert len(result) == 1
    assert all(len(r['ubiquitin_history']) == 2 for r in result)
    assert all(len(r['reaction_history']) == 2 for r in result)
    assert all(len(r['donor_history']) == 2 for r in result)
    assert all(len(r['context_history']) == 2 for r in result)

def test_simulate_deprot_steps_output_structure_nested_with_SMAC():
    """
    Test that simulate_deprot_steps returns correct structure on nested dictionary. 
    
    If there is a SMAC to deprotect so len(result) should be 2
    """
    result = simulate_deprot_steps(history_dict_3_3a)

    assert len(result) == 2
    assert all(len(r['ubiquitin_history']) == 7 for r in result)
    assert all(len(r['reaction_history']) == 7 for r in result)
    assert all(len(r['donor_history']) == 7 for r in result)
    assert all(len(r['context_history']) == 7 for r in result)

def test_simulate_deprot_steps_output_structure_nested_without_SMAC():
    """
    Test that simulate_deprot_steps returns correct structure on nested dictionary. 
    
    If there is no SMAC to deprotect so len(result) should be 1
    """
    result = simulate_deprot_steps(history_dict_3_3b)

    assert len(result) == 1
    assert all(len(r['ubiquitin_history']) == 8 for r in result)
    assert all(len(r['reaction_history']) == 8 for r in result)
    assert all(len(r['donor_history']) == 8 for r in result)
    assert all(len(r['context_history']) == 8 for r in result)

def test_simulate_deprot_steps_appends_correct_reactions():
    """Check that FAKE_deprot is only appended correctly as there is no SMAC on the species."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1_K48_aboc)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    result = simulate_deprot_steps(history_dict)
    assert [r['reaction_history'][-1] for r in result] == ['SMAC_deprot', 'FAKE_deprot']

def test_simulate_deprot_steps_appends_correct_reactions():
    """Check that SMAC_deprot and FAKE_deprot are appended correctly."""
    result = simulate_deprot_steps(history_dict_3_3a)
    assert [r['reaction_history'][-1] for r in result] == ['SMAC_deprot', 'FAKE_deprot']

def test_simulate_deprot_steps_multimer_contains_expected_keys():
    """Verify new multimers contain expected structure."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1_K63_aboc)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    result = simulate_deprot_steps(history_dict)
    for r in result:
        assert 'protein' in r['ubiquitin_history'][-1]
        assert 'branching_sites' in r['ubiquitin_history'][-1]

def test_simulate_deprot_steps_context_structure():
    """Ensure context entries include necessary fields."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    result = simulate_deprot_steps(history_dict)
    for r in result:
        assert 'chain_number_list' in r['context_history'][-1]
        assert 'chain_length_list' in r['context_history'][-1]

def test_simulate_deprot_steps_missing_keys_raises_error():
    """Check missing keys raise KeyError."""
    history_dict = {
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': ['']
    }
    with pytest.raises(KeyError):
        simulate_deprot_steps(history_dict)

# =====================================================================
# =================== Combined E2 + Deprot Step Test ==================
# =====================================================================

def test_simulate_e2_then_deprot_combination():
    """Test combined simulation: E2 step followed by deprot step for each result."""
    acceptor, context = iterate_through_ubiquitin(histag_ubi_ubq_1)
    history_dict = {
        'ubiquitin_history': [acceptor],
        'reaction_history': [''],
        'donor_history': [''],
        'context_history': [context]
    }
    donor_list = [ubi_ubq_1_K48_SMAC, ubi_ubq_1_K63_SMAC]
    e2_results = simulate_E2_steps(history_dict, donor_list)

    all_deprot_results = []
    for intermediate in e2_results:
        deprot_results = simulate_deprot_steps(intermediate)
        all_deprot_results.extend(deprot_results)

    # Expect 2 E2 results (only 2 reactions allowed) x 2 deprot each = 4 total
    assert len(all_deprot_results) == 4
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
    assert result == "Ubc13/Mms2 (branching)"


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
        [1, "K48", 2], [2, "K63", 3], [3, "K48", 4]
    ]
    new_bound_lysine = [3, "K48", 4]
    product_conjugated_lysines.append([new_bound_lysine])
    result = determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert result == "elongation"


def test_determine_elongation_or_branching_branching():
    """Test when new bound lysine appears twice: should return 'branching'."""
    product_conjugated_lysines = [
        [1, "K48", 2], [2, "K48", 3], [1, "K63", 3]
    ]
    new_bound_lysine = [1, "K63", 3]
    product_conjugated_lysines.append([new_bound_lysine])
    result = determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert result == "branching"


def test_determine_elongation_or_branching_invalid_count_zero():
    """Test when new bound lysine does not exist at all: should raise TypeError."""
    product_conjugated_lysines = [
        [1, "K48", 2], [2, "K63", 3]
    ]
    new_bound_lysine = [3, "K48", 4]
    product_conjugated_lysines.append([new_bound_lysine])
    with pytest.raises(TypeError) as exc_info:
        determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert "expected 1 or 2" in str(exc_info.value)


def test_determine_elongation_or_branching_invalid_count_more_than_two():
    """Test when new bound lysine appears more than twice: should raise TypeError."""
    product_conjugated_lysines = [
        [1, "K63", 2], [2, "K48", 3], [2, "K63", 4], [2, "K29", 5]
    ]
    new_bound_lysine = [2, "K63", 4]
    product_conjugated_lysines.append([new_bound_lysine])
    with pytest.raises(TypeError) as exc_info:
        determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert "expected 1 or 2" in str(exc_info.value)


def test_determine_elongation_or_branching_with_empty_product_list():
    """Test when product_conjugated_lysines is length 1 and new lysine is added once."""
    product_conjugated_lysines = [[1, "K48", 2]]
    new_bound_lysine = [1, "K48", 2]
    product_conjugated_lysines.append([new_bound_lysine])
    result = determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine)
    assert result == "elongation"

# =====================================
# Tests for determine_reaction_type
# =====================================

def test_determine_reaction_type_k48():
    """Test that K48 lysine returns correct reaction type."""
    new_bound_lysine = [3, "K48", 4]
    result = determine_reaction_type(new_bound_lysine)
    assert result == "K48_reaction"


def test_determine_reaction_type_k63():
    """Test that K63 lysine returns correct reaction type."""
    new_bound_lysine = [5, "K63", 6]
    result = determine_reaction_type(new_bound_lysine)
    assert result == "K63_reaction"


def test_determine_reaction_type_invalid_lysine_raises():
    """Test that an invalid lysine site raises TypeError."""
    new_bound_lysine = [2, "K33", 3]
    with pytest.raises(TypeError) as exc_info:
        determine_reaction_type(new_bound_lysine)
    assert "does not contain K48 or K63" in str(exc_info.value)


def test_determine_reaction_type_non_string_site_raises():
    """Test that a non-string lysine site raises TypeError."""
    new_bound_lysine = [4, 48, 5]
    with pytest.raises(TypeError) as exc_info:
        determine_reaction_type(new_bound_lysine)
    assert "does not contain K48 or K63" in str(exc_info.value)

# =====================================
# Tests for determine_Ube2K_elongation
# =====================================

@pytest.mark.parametrize(
    "product_conjugated_lysines, new_bound_lysine, expected_result, expected_exception, note",
    [
        # Test 1: K63 linkage → should return Ube2K
        (
            [[1, "K48", 2], [2, "K63", 3], [3, "K48", 4]],
            [3, "K48", 4],
            "Ube2K",
            None,
            "Test when previous linkage was K63 → should return Ube2K."
        ),

        # Test 2: K48 linkage → should return gp78/Ube2g2
        (
            [[1, "K48", 2], [2, "K48", 3], [3, "K48", 4]],
            [3, "K48", 4],
            "gp78/Ube2g2",
            None,
            "Test when previous linkage was K48 → should return gp78/Ube2g2."
        ),

        # Test 4: Multiple linkages present, correct K63 found → Ube2K
        (
            [[1, "K48", 2], [2, "K48", 3], [3, "K63", 4], [4, "K48", 5], [2, "K63", 6]],
            [4, "K48", 5],
            "Ube2K",
            None,
            "Test when multiple linkages exist but correct K63 linkage is found → should return Ube2K."
        ),

        # Test 5: Multiple linkages present, correct K48 found → gp78/Ube2g2
        (
            [[1, "K48", 2], [2, "K48", 3], [3, "K48", 4], [4, "K48", 5], [3, "K63", 6]],
            [4, "K48", 5],
            "gp78/Ube2g2",
            None,
            "Test when multiple linkages exist but correct K48 linkage is found → should return gp78/Ube2g2."
        ),
    ]
)
def test_determine_Ube2K_elongation_parametrized(
    product_conjugated_lysines,
    new_bound_lysine,
    expected_result,
    expected_exception,
    note
):
    """Parametrized test for determine_Ube2K_elongation — covers K63, K48, missing, multiple linkages."""
    if expected_exception:
        with pytest.raises(expected_exception) as exc_info:
            determine_Ube2K_elongation(product_conjugated_lysines, new_bound_lysine)
        # Optional: can assert error message if desired
    else:
        result = determine_Ube2K_elongation(product_conjugated_lysines, new_bound_lysine)
        assert result == expected_result

# =====================================
# Tests for validate_conjugated_lysines
# =====================================

def test_validate_conjugated_lysines_all_valid():
    """Test that validation passes when only K48 and K63 are present."""
    context = {
        "conjugated_lysines": [[1, "K48", 2], [2, "K63", 3], [3, "K48", 4]]
    }
    # Should not raise any exceptions
    validate_conjugated_lysines(context)


def test_validate_conjugated_lysines_with_invalid_entry():
    """Test that validation raises TypeError when unsupported lysines are included."""
    context = {
        "conjugated_lysines": [[1, "K48", 2], [2, "K33", 3], [3, "K63", 4]]
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

def test_validate_conjugated_lysines_raises_on_malformed_entries():
    """
    Test that validate_conjugated_lysines raises TypeError when conjugated_lysines
    contains only malformed entries (not lists of length 3).
    """
    context = {
        "conjugated_lysines": [[1, "K48"], [2, "K48", 3], [3, "K63", 4]]
    }

    with pytest.raises(TypeError) as exc_info:
        validate_conjugated_lysines(context)

    assert "malformed conjugated_lysines" in str(exc_info.value)


# ============================================================
# Build reaction sequence for testing assign_correct_E2_enzyme
# ============================================================
# ===== Histories are there in-case you want to test later ===
# === First Reaction Path ===

acceptor_1_0, context_1_0 = iterate_through_ubiquitin(histag_ubi_ubq_1_K63_aboc)
history_dict_1_0 = {}
history_dict_1_0['ubiquitin_history'] = [acceptor_1_0]
history_dict_1_0['reaction_history'] = ['']
history_dict_1_0['donor_history'] = ['']
history_dict_1_0['context_history'] = [context_1_0]

# Dimer formation (K48 elongation)
acceptor = acceptor_1_0.copy()
context = context_1_0.copy()
donor = ubi_ubq_1_K48_SMAC
reaction = "K48"
acceptor_1_1a, context_1_1a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_1_1a = {}
history_dict_1_1a['ubiquitin_history'] = history_dict_1_0['ubiquitin_history'] + [acceptor_1_1a]
history_dict_1_1a['reaction_history'] = history_dict_1_0['reaction_history'] + [reaction]
history_dict_1_1a['donor_history'] = history_dict_1_0['donor_history'] + [donor]
history_dict_1_1a['context_history'] = history_dict_1_0['context_history'] + [context_1_1a]

# Dimer deprotection
acceptor, context = acceptor_1_1a.copy(), context_1_1a.copy()
acceptor_1_1b, context_1_1b = ubiquitin_simulation(acceptor, '', "FAKE_deprot")
history_dict_1_1b = {}
history_dict_1_1b['ubiquitin_history'] = history_dict_1_1a['ubiquitin_history'] + [acceptor_1_1b]
history_dict_1_1b['reaction_history'] = history_dict_1_1a['reaction_history'] + [reaction]
history_dict_1_1b['donor_history'] = history_dict_1_1a['donor_history'] + [donor]
history_dict_1_1b['context_history'] = history_dict_1_1a['context_history'] + [context_1_1b]


# Trimer formation (K63 elongation)
acceptor, context = acceptor_1_1b.copy(), context_1_1b.copy()
donor = ubi_ubq_1_K48_ABOC_K63_ABOC
reaction = "K63"
acceptor_1_2a, context_1_2a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_1_2a = {}
history_dict_1_2a['ubiquitin_history'] = history_dict_1_1b['ubiquitin_history'] + [acceptor_1_2a]
history_dict_1_2a['reaction_history'] = history_dict_1_1b['reaction_history'] + [reaction]
history_dict_1_2a['donor_history'] = history_dict_1_1b['donor_history'] + [donor]
history_dict_1_2a['context_history'] = history_dict_1_1b['context_history'] + [context_1_2a]


# Trimer deprotection
acceptor, context = acceptor_1_2a.copy(), context_1_2a.copy()
acceptor_1_2b, context_1_2b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")
history_dict_1_2b = {}
history_dict_1_2b['ubiquitin_history'] = history_dict_1_2a['ubiquitin_history'] + [acceptor_1_2b]
history_dict_1_2b['reaction_history'] = history_dict_1_2a['reaction_history'] + [reaction]
history_dict_1_2b['donor_history'] = history_dict_1_2a['donor_history'] + [donor]
history_dict_1_2b['context_history'] = history_dict_1_2a['context_history'] + [context_1_2b]

# Tetramer formation (K48 branching)
acceptor, context = acceptor_1_2b.copy(), context_1_2b.copy()
donor = ubi_ubq_1_K48_ABOC_K63_ABOC
reaction = "K48"
acceptor_1_3a, context_1_3a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_1_3a = {}
history_dict_1_3a['ubiquitin_history'] = history_dict_1_2b['ubiquitin_history'] + [acceptor_1_3a]
history_dict_1_3a['reaction_history'] = history_dict_1_2b['reaction_history'] + [reaction]
history_dict_1_3a['donor_history'] = history_dict_1_2b['donor_history'] + [donor]
history_dict_1_3a['context_history'] = history_dict_1_2b['context_history'] + [context_1_3a]

# Tetramer deprotection
acceptor, context = acceptor_1_3a.copy(), context_1_3a.copy()
acceptor_1_3b, context_1_3b = ubiquitin_simulation(acceptor, '', "FAKE_deprot")
history_dict_1_3b = {}
history_dict_1_3b['ubiquitin_history'] = history_dict_1_3a['ubiquitin_history'] + [acceptor_1_3b]
history_dict_1_3b['reaction_history'] = history_dict_1_3a['reaction_history'] + [reaction]
history_dict_1_3b['donor_history'] = history_dict_1_3a['donor_history'] + [donor]
history_dict_1_3b['context_history'] = history_dict_1_3a['context_history'] + [context_1_3b]



# === Second Reaction Path ===

acceptor_2_0, context_2_0 = iterate_through_ubiquitin(histag_ubi_ubq_1_K63_aboc)
history_dict_2_0 = {}
history_dict_2_0['ubiquitin_history'] =  [acceptor_2_0]
history_dict_2_0['reaction_history'] =  ['']
history_dict_2_0['donor_history'] =  ['']
history_dict_2_0['context_history'] =  [context_2_0]

# Dimer formation (K48 elongation)
acceptor, context = acceptor_2_0.copy(), context_2_0.copy()
donor = ubi_ubq_1_K48_ABOC_K63_SMAC
reaction = "K48"
acceptor_2_1a, context_2_1a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_2_1a = {}
history_dict_2_1a['ubiquitin_history'] = history_dict_2_0['ubiquitin_history'] + [acceptor_2_1a]
history_dict_2_1a['reaction_history'] = history_dict_2_0['reaction_history'] + [reaction]
history_dict_2_1a['donor_history'] = history_dict_2_0['donor_history'] + [donor]
history_dict_2_1a['context_history'] = history_dict_2_0['context_history'] + [context_2_1a]

# Dimer deprotection
acceptor, context = acceptor_2_1a.copy(), context_2_1a.copy()
acceptor_2_1b, context_2_1b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")
history_dict_2_1b = {}
history_dict_2_1b['ubiquitin_history'] = history_dict_2_1a['ubiquitin_history'] + [acceptor_2_1b]
history_dict_2_1b['reaction_history'] = history_dict_2_1a['reaction_history'] + [reaction]
history_dict_2_1b['donor_history'] = history_dict_2_1a['donor_history'] + [donor]
history_dict_2_1b['context_history'] = history_dict_2_1a['context_history'] + [context_2_1b]


# Trimer formation (K63 elongation)
acceptor, context = acceptor_2_1b.copy(), context_2_1b.copy()
donor = ubi_ubq_1_K63_SMAC
reaction = "K63"
acceptor_2_2a, context_2_2a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_2_2a = {}
history_dict_2_2a['ubiquitin_history'] = history_dict_2_1b['ubiquitin_history'] + [acceptor_2_2a]
history_dict_2_2a['reaction_history'] = history_dict_2_1b['reaction_history'] + [reaction]
history_dict_2_2a['donor_history'] = history_dict_2_1b['donor_history'] + [donor]
history_dict_2_2a['context_history'] = history_dict_2_1b['context_history'] + [context_2_2a]

# Trimer deprotection
acceptor, context = acceptor_2_2a.copy(), context_2_2a.copy()
acceptor_2_2b, context_2_2b = ubiquitin_simulation(acceptor, '', "FAKE_deprot")
history_dict_2_2b = {}
history_dict_2_2b['ubiquitin_history'] = history_dict_2_2a['ubiquitin_history'] + [acceptor_2_2b]
history_dict_2_2b['reaction_history'] = history_dict_2_2a['reaction_history'] + [reaction]
history_dict_2_2b['donor_history'] = history_dict_2_2a['donor_history'] + [donor]
history_dict_2_2b['context_history'] = history_dict_2_2a['context_history'] + [context_2_2b]

# Tetramer formation (K48 Ube2K elongation as previous is K63)
acceptor, context = acceptor_2_2b.copy(), context_2_2b.copy()
donor = ubi_ubq_1_K48_ABOC_K63_ABOC
reaction = "K48"
acceptor_2_3a, context_2_3a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_2_3a = {}
history_dict_2_3a['ubiquitin_history'] = history_dict_2_2b['ubiquitin_history'] + [acceptor_2_3a]
history_dict_2_3a['reaction_history'] = history_dict_2_2b['reaction_history'] + [reaction]
history_dict_2_3a['donor_history'] = history_dict_2_2b['donor_history'] + [donor]
history_dict_2_3a['context_history'] = history_dict_2_2b['context_history'] + [context_2_3a]

# Tetramer deprotection
acceptor, context = acceptor_2_3a.copy(), context_2_3a.copy()
acceptor_2_3b, context_2_3b = ubiquitin_simulation(acceptor, '', "SMAC_deprot")
history_dict_2_3b = {}
history_dict_2_3b['ubiquitin_history'] = history_dict_2_3a['ubiquitin_history'] + [acceptor_2_3b]
history_dict_2_3b['reaction_history'] = history_dict_2_3a['reaction_history'] + [reaction]
history_dict_2_3b['donor_history'] = history_dict_2_3a['donor_history'] + [donor]
history_dict_2_3b['context_history'] = history_dict_2_3a['context_history'] + [context_2_3b]

# Pentamer formation (K63 branching)
acceptor, context = acceptor_2_3b.copy(), context_2_3b.copy()
donor = ubi_ubq_1_K48_ABOC_K63_ABOC
reaction = "K63"
acceptor_2_4a, context_2_4a = ubiquitin_simulation(acceptor, donor, reaction)
history_dict_2_4a = {}
history_dict_2_4a['ubiquitin_history'] = history_dict_2_3b['ubiquitin_history'] + [acceptor_2_4a]
history_dict_2_4a['reaction_history'] = history_dict_2_3b['reaction_history'] + [reaction]
history_dict_2_4a['donor_history'] = history_dict_2_3b['donor_history'] + [donor]
history_dict_2_4a['context_history'] = history_dict_2_3b['context_history'] + [context_2_4a]

# Pentamer deprotection
acceptor, context = acceptor_2_4a.copy(), context_2_4a.copy()
acceptor_2_4b, context_2_4b = ubiquitin_simulation(acceptor, '', "FAKE_deprot")
history_dict_2_4b = {}
history_dict_2_4b['ubiquitin_history'] = history_dict_2_4a['ubiquitin_history'] + [acceptor_2_4b]
history_dict_2_4b['reaction_history'] = history_dict_2_4a['reaction_history'] + [reaction]
history_dict_2_4b['donor_history'] = history_dict_2_4a['donor_history'] + [donor]
history_dict_2_4b['context_history'] = history_dict_2_4a['context_history'] + [context_2_4b]

# Parameterized test for assign_correct_E2_enzyme
@pytest.mark.parametrize(
    "reactant, product, expected_enzyme, note",
    [
        (context_1_0, context_1_1a, "gp78/Ube2g2", "Dimer formation via K48 results in elongation"),
        (context_1_1b, context_1_2a, "Ubc13/Mms2", "Trimer formation via K63 results in elongation"),
        (context_1_2b, context_1_3a, "Ube2K", "Tetramer formation via K48 results in branching"),
        (context_2_0, context_2_1a, "gp78/Ube2g2", "Dimer formation via K48 results in elongation (2nd path)"),
        (context_2_1b, context_2_2a, "Ubc13/Mms2", "Trimer formation via K63 results in elongation (2nd path)"),
        (context_2_2b, context_2_3a, "Ube2K", "Tetramer formation via K48 results in elongation but (Ube2K elongation) (2nd path)"),
        (context_2_3b, context_2_4a, "Ubc13/Mms2 (branching)", "Pentamer formation via K63 results in branching"),
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
        "max_chain_number": 2,
        "conjugated_lysines": [[1, "K48", 2]],
        "chain_number_list": [1],
        "chain_length_list": []
    }
    with pytest.raises(TypeError) as exc_info:
        assign_correct_E2_enzyme(context, context)
    assert "product_max_chain_number: 2 != reactant_max_chain_number + 1: 3" in str(exc_info.value)