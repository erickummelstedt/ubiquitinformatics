import types
import pytest
import json
import copy
from copy import deepcopy
import sys
import logging
from typing import Dict, List, Any, Union
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
    inject_branching_sites, \
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
    

def simulate_E2_step(
    history_dict: dict,
    donor_list: list
) -> list[dict]:
    """
    Simulate reactions for a donor list at K48 or K63 sites.

    Args:
        history_dict (dict): Dictionary containing:
            - 'ubiquitin_history' (list): List of previously accepted protein states.
            - 'reaction_history' (list): List of past reactions applied.
            - 'donor_history' (list): List of ubiquitin monomers used.
            - 'context_history' (list): List of contexts.
        donor_list (list): Available monomers to test in reactions.

    Returns:
        list[dict]: List of new history dictionaries generated from the reactions.
    """
    reaction_types = ['K48', 'K63']
    results = []

    for reaction in reaction_types:
        for monomer in donor_list:
            # Skip incompatible donor-reaction pairs
            # This is not allowed because the donor has a protecting group at the same site it would conjugate through.
            # Allowing this would result in a self-reaction (monomer reacting with itself), which is invalid.
            # if the K48 has not protecting group it will react itself in a K48 reaction.
            if (
                ((monomer == ubi_ubq_1_K48_SMAC) and (reaction == "K63")) or
                ((monomer == ubi_ubq_1_K63_SMAC) and (reaction == "K48"))
            ):
                continue
            last_acceptor = history_dict['ubiquitin_history'][-1]
            new_multimer, new_context = ubiquitin_simulation(last_acceptor, monomer, reaction)

            results.append({
                'ubiquitin_history': history_dict['ubiquitin_history'] + [new_multimer],
                'reaction_history': history_dict['reaction_history'] + [reaction],
                'donor_history': history_dict['donor_history'] + [monomer],
                'context_history': history_dict['context_history'] + [new_context]
            })

    return results


def simulate_deprot_step(history_dict: dict) -> list[dict]:
    """
    Simulate SMAC deprotection and buffer wash reactions.

    Args:
        history_dict (dict): Dictionary containing:
            - 'ubiquitin_history' (list): List of previously accepted protein states.
            - 'reaction_history' (list): List of past reactions applied.
            - 'donor_history' (list): List of ubiquitin monomers used.
            - 'context_history' (list): List of contexts.

    Returns:
        list[dict]: List of new history dictionaries generated from the reactions.
    """
    reaction_list = ['SMAC_deprot', 'FAKE_deprot']
    results = []

    for reaction in reaction_list:
        last_acceptor = history_dict['ubiquitin_history'][-1]
        new_multimer, new_context = ubiquitin_simulation(last_acceptor, '', reaction)

        results.append({
            'ubiquitin_history': history_dict['ubiquitin_history'] + [new_multimer],
            'reaction_history': history_dict['reaction_history'] + [reaction],
            'donor_history': history_dict['donor_history'] + [''],
            'context_history': history_dict['context_history'] + [new_context]
        })

    return results


def assign_enzyme(reaction, elongation_or_branching):
    """
    Assigns the correct E2 enzyme based on reaction type and linkage pattern.

    Args:
        reaction (str): Type of reaction ('K48_reaction' or 'K63_reaction').
        elongation_or_branching (str): 'elongation' or 'branching'.

    Returns:
        str: Name of the correct E2 enzyme.

    Raises:
        TypeError: If combination is invalid.
    """
    if reaction == "K48_reaction" and elongation_or_branching == "elongation":
        return "gp78/Ube2g2"
    elif reaction == "K48_reaction" and elongation_or_branching == "branching":
        return "Ube2K"
    elif reaction == "K63_reaction":
        return "Ubc13/Mms2"
    else:
        raise TypeError(
            f"Invalid combination of reaction: {reaction} and "
            f"elongation_or_branching: {elongation_or_branching}"
        )

def determine_elongation_or_branching(product_conjugated_lysines, new_bound_lysine):
    """
    Determines whether the addition of a ubiquitin results in elongation or branching
    based on how many times the same chain number appears.

    Args:
        product_conjugated_lysines (list): List of all conjugated lysines in the product context.
        new_bound_lysine (list): The site where the most recent ubiquitin was added (e.g. [chain_number, lysine_site]).

    Returns:
        str: 'elongation' if the chain appears once, 'branching' if it appears twice.

    Raises:
        TypeError: If the count of the chain number is not 1 or 2.
    """
    target_chain = new_bound_lysine[0]

    count = sum(1 for entry in product_conjugated_lysines if entry[0] == target_chain)

    if count == 1:
        return "elongation"
    elif count == 2:
        return "branching"
    else:
        raise TypeError(
            f"Count of chain number '{target_chain}' is {count}, expected 1 or 2"
        )

def determine_reaction_type(new_bound_lysine):
    """
    Determines the type of reaction (K48 or K63) from the lysine site of the newly bound ubiquitin.

    Args:
        new_bound_lysine (list): A list with [chain_number, lysine_site], e.g., [3, "K48"]

    Returns:
        str: 'K48_reaction' or 'K63_reaction'

    Raises:
        TypeError: If the lysine site is not K48 or K63.
    """
    lysine_site = new_bound_lysine[1]

    if lysine_site == "K48":
        return "K48_reaction"
    elif lysine_site == "K63":
        return "K63_reaction"
    else:
        raise TypeError(
            f"new_bound_lysine: {new_bound_lysine} "
            f"does not contain K48 or K63"
        )

def validate_conjugated_lysines(context):
    """
    Validates that the conjugated lysines in the context are only K48 or K63.

    Args:
        context (dict): A context dictionary containing a 'conjugated_lysines' list.

    Raises:
        TypeError: If unsupported lysines are found.
    """
    lysine_types = {
        lys[1] for lys in context["conjugated_lysines"]
        if isinstance(lys, list) and len(lys) == 2
    }
    unsupported_lysines = lysine_types - {"K48", "K63"}
    if unsupported_lysines:
        raise TypeError(
            f"Context contains unsupported conjugated lysines: {unsupported_lysines}. "
            "This function only supports K48 and K63 conjugated lysines."
        )
    
def assign_correct_E2_enzyme(
        reactant_context: str | dict,
        product_context: str | dict
        ):
    """
    Determine the enzyme used in a reaction based on structural context changes
    between a reactant and product ubiquitin chain.

    Args:
        reactant_context (str | dict): Reactant context dictionary or JSON string.
        product_context (str | dict): Product context dictionary or JSON string.

    Returns:
        str: Enzyme type ('Ube2K', 'Ube2g2', or 'Ubc13/Mms2').
    """
    reactant_context = convert_json_to_dict(reactant_context)
    product_context = convert_json_to_dict(product_context)

    if int(product_context["max_chain_number"]) != int(reactant_context["max_chain_number"]) + 1:
        raise TypeError(
            f"product_max_chain_number: {product_context['max_chain_number']} != "
            f"reactant_max_chain_number + 1: {int(reactant_context['max_chain_number']) + 1}"
        )

    validate_conjugated_lysines(reactant_context)
    validate_conjugated_lysines(product_context)

    reactant_conjugated_lysines = reactant_context["conjugated_lysines"].copy() + [[]]
    product_conjugated_lysines = product_context["conjugated_lysines"].copy()

    for _, (reactant_lysine, product_lysine) in enumerate(
        zip(reactant_conjugated_lysines, product_conjugated_lysines)
    ):
        if reactant_lysine != product_lysine:
            new_bound_lysine = product_lysine
            break
    else:
        raise TypeError("No new conjugation site detected between reactant and product contexts.")

    reaction = determine_reaction_type(new_bound_lysine)

    elongation_or_branching = determine_elongation_or_branching(
        product_conjugated_lysines, new_bound_lysine
    )

    enzyme = assign_enzyme(reaction, elongation_or_branching)

    return enzyme