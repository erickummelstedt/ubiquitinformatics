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
    find_max_chain_number, \
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




def find_Ube2K_or_Uce2g2(
        reactant_context: str|dict,
        product_context: str|dict
        ):
    """
    Find the enzyme type based on the reactant and product dictionaries.
    Args:
        reactant_dictionary (str|dict): Reactant dictionary or JSON string.
        product_dictionary (str|dict): Product dictionary or JSON string.
    Returns:
        str: Enzyme type ('Ube2K' or 'Ube2g2').
    """

    # Convert JSON strings to dictionaries if necessary
    reactant_context = convert_json_to_dict(reactant_context)    
    product_context = convert_json_to_dict(product_context)

    reactant_chain_number_list  = reactant_context["chain_number_list"].copy()
    reactant_chain_length_list = reactant_context["chain_length_list"].copy()
    reactant_free_lysine_list = reactant_context["free_lysines"].copy()
    reactant_bound_lysine_list = reactant_context["conjugated_lysines"].copy()
    reactant_max_chain_number = reactant_context["max_chain_number"]

    product_chain_number_list = product_context["chain_number_list"].copy()
    product_chain_length_list = product_context["chain_length_list"].copy()
    product_free_lysine_list = product_context["free_lysines"].copy()
    product_bound_lysine_list = product_context["conjugated_lysines"].copy()
    product_max_chain_number = product_context["max_chain_number"]

    final_dictionary = {'reactant_chain_numbers' : reactant_chain_number_list, 
                        'reactant_chain_lengths' : reactant_chain_length_list, 
                        'reactant_free_lysines' : reactant_free_lysine_list, 
                        'reactant_bound_lysines' : reactant_bound_lysine_list, 
                        'reactant_max_chain_number' : reactant_max_chain_number, 
                        'product_chain_numbers' : product_chain_number_list,
                        'product_chain_lengths' : product_chain_length_list,
                        'product_free_lysines' : product_free_lysine_list,
                        'product_bound_lysines' : product_bound_lysine_list,
                        'product_max_chain_number' : product_max_chain_number
    }

    reactant_bound_lysines = final_dictionary['reactant_bound_lysines'].copy()
    reactant_bound_lysines = reactant_bound_lysines + [[]]

    reactant_max_chain_number = int(final_dictionary['reactant_max_chain_number'])

    product_bound_lysines = final_dictionary['product_bound_lysines'].copy()

    ## maybe try and break...
    for index, (x,y) in enumerate(zip(reactant_bound_lysines,product_bound_lysines)):
        ## if it is ie a monomer and we are acceptorising on K48 it is not branching
        ## add in max chain number is == 1 then it is a acceptorisation and hence gp78/Ube2g2
        if  (x!=y) & (index == 0) & (reactant_max_chain_number == 1):
            last_bound_site__before_new_ubi = [1,'K48']
            break
        ## add in max chain number is != 1 then it is not a acceptorsation and hence Ube2K
        elif (x!=y) & (index == 0) & (reactant_max_chain_number != 1):
            last_bound_site__before_new_ubi = [1,'K63']
            break
        elif x==y:
            #print(index, (x,y))
            last_bound_site__before_new_ubi = y
        else:
            break
    
    if last_bound_site__before_new_ubi[1] == 'K48':
        enzyme = 'gp78/Ube2g2'
    else: 
        enzyme = 'Ube2K'
    return enzyme









def find_K63_branching(reactant_dictionary, product_dictionary):
    global chain_number_list
    global chain_length_list
    global free_lysine_list
    global bound_lysine_list
    
    ## do everything with reactant and product...
    chain_number_list = [1]
    chain_length_list = []
    free_lysine_list = []
    bound_lysine_list = []
    
    if isinstance(reactant_dictionary, str):
        x = reactant_dictionary.replace("'", "\"")
        reactant_dictionary = json.loads(x)
    
    if isinstance(product_dictionary, str):
        x = product_dictionary.replace("'", "\"")
        product_dictionary = json.loads(x)

    logging.info(chain_number_list)
    adapted_dictionary = inner_wrapper_find_K63_branching(reactant_dictionary)

    reactant_chain_number_list  = chain_number_list.copy()
    reactant_chain_length_list = chain_length_list.copy()
    reactant_free_lysine_list = free_lysine_list.copy()
    reactant_bound_lysine_list = bound_lysine_list.copy()
    reactant_max_chain_number = chain_number_list[-1]-1

    chain_number_list = [1]
    chain_length_list = []
    free_lysine_list = []
    bound_lysine_list = []


    logging.info(chain_number_list)
    adapted_dictionary = inner_wrapper_find_K63_branching(product_dictionary)

    product_chain_number_list = chain_number_list.copy()
    product_chain_length_list = chain_length_list.copy()
    product_free_lysine_list = free_lysine_list.copy()
    product_bound_lysine_list = bound_lysine_list.copy()
    product_max_chain_number = chain_number_list[-1]-1



    final_dictionary = {'reactant_chain_numbers' : reactant_chain_number_list, 
                        'reactant_chain_lengths' : reactant_chain_length_list, 
                        'reactant_free_lysines' : reactant_free_lysine_list, 
                        'reactant_bound_lysines' : reactant_bound_lysine_list, 
                        'reactant_max_chain_number' : reactant_max_chain_number, 
                        'product_chain_numbers' : product_chain_number_list,
                        'product_chain_lengths' : product_chain_length_list,
                        'product_free_lysines' : product_free_lysine_list,
                        'product_bound_lysines' : product_bound_lysine_list,
                        'product_max_chain_number' : product_max_chain_number
    }

    reactant_bound_lysines = final_dictionary['reactant_bound_lysines'].copy()
    reactant_bound_lysines = reactant_bound_lysines + [[]]

    reactant_max_chain_number = int(final_dictionary['reactant_max_chain_number'])

    product_bound_lysines = final_dictionary['product_bound_lysines'].copy()

    ## maybe try and break...
    for index, (x,y) in enumerate(zip(reactant_bound_lysines,product_bound_lysines)):
        ## if it is ie a monomer and we are acceptorising on K48 it is not branching
        ## add in max chain number is == 1 then it is a acceptorisation and hence gp78/Ube2g2
        if  (x!=y) & (index == 0) & (reactant_max_chain_number == 1):
            last_bound_site__before_new_ubi = [1,'K63']
            break
        ## add in max chain number is != 1 then it is not a acceptorsation and hence Ube2K
        elif (x!=y) & (index == 0) & (reactant_max_chain_number != 1):
            last_bound_site__before_new_ubi = [1,'K48']
            break
        elif x==y:
            #print(index, (x,y))
            last_bound_site__before_new_ubi = y
        else:
            break
    
    if last_bound_site__before_new_ubi[1] == 'K63':
        enzyme = 'Ube13/Mms2'
    else: 
        enzyme = 'Ube13/Mms2_branching'
    return enzyme

    #if not(specific_ubi_num in chain_number_list[:-1]): $
    #    raise TypeError('the ubiquitin specific is outside the range of the ubiquitin multimer')
