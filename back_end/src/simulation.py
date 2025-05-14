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
    

def simulate_E2_steps(
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
            # Skip incompatible donor-reaction pairs; the self-reaction is not allowed.
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

            # If one the chain number does not increase by 1, then no ubiquitin was added or more than 1 ubiquitin was added 
            last_context = history_dict['context_history'][-1]
            if int(new_context['max_chain_number']) != (int(last_context['max_chain_number']) + 1):
                continue

            # LAST TODO 
            # Apply the enzyme choice using new context and old context

            results.append({
                'ubiquitin_history': history_dict['ubiquitin_history'] + [new_multimer],
                'reaction_history': history_dict['reaction_history'] + [reaction],
                'donor_history': history_dict['donor_history'] + [monomer],
                'context_history': history_dict['context_history'] + [new_context]
            })

    return results


def simulate_deprot_steps(history_dict: dict) -> list[dict]:
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

        # Check if the reaction is SMAC and if there is no change in the reaction
        if reaction == 'SMAC_deprot':
            # Check if the reaction is SMAC and if there is no change in the reaction
            # If so, skip this reaction
            if str(new_multimer) == str(last_acceptor):
                continue
        
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


# create synthesis database
# input are only; acceptors and donors 
# always start with E2 reaction
# output should be compared to globally deprotected acceptors

## instead of list of lists use a list of dictionaries
def create_synthesis_database(
        acceptors : list, 
        donors : list
        ):
    """
    Create synthesis dataframes for acceptors and donors.
    Args:
        acceptors (list): List of acceptor proteins.
        donors (list): List of donor proteins.          
    Returns:
        tuple: Dataframes for acceptors, reactions, and monomers.
    """

    # create history dictionaries
    list_of_history_dicts = []
    for i in acceptors:
        acceptor, context = iterate_through_ubiquitin(i)
        history_dict = {
            'ubiquitin_history': [acceptor],
            'reaction_history': [''],
            'donor_history': [''],
            'context_history': [context]
        }
        list_of_history_dicts.append(history_dict)


    # simulate E2 steps on histories
    for i in list_of_history_dicts:
        list_of_history_dicts = simulate_E2_steps(i, donors)

    # simulate deprot steps on histories
    for i in list_of_history_dicts:
        list_of_history_dicts = simulate_deprot_steps(i, donors)
        



    
    
    
    
    ## starting lists
    #growing_acceptor_history_list = [[ubi_ubq_1_K48_SMAC], [ubi_ubq_1_K63_SMAC], [ubi_ubq_1_K48_SMAC_K63_ABOC], [ubi_ubq_1_K48_ABOC_K63_SMAC], [ubi_ubq_1_K48_ABOC_K63_ABOC]]
    ## take list from dataframe... 
    growing_acceptor_history_list = input_acceptor_history_list.copy()
    growing_reaction_history_list = input_reaction_history_list.copy()
    growing_monomer_history_list = input_monomer_history_list.copy()

    new_growing_acceptor_history_list = []
    new_growing_reaction_history_list = []
    new_growing_monomer_history_list = []

    if enzyme_of_deprot_first == 'enzyme': 
        #for i in range(5):
        #for i in range(3):
        for (a,b,c) in zip(growing_acceptor_history_list, growing_reaction_history_list, growing_monomer_history_list):
            ## get the next acceptor...
            current_acceptor_history, current_reaction_history, current_monomer_history = simulate_reactions_step(a,b,c, ubi_donor_list)
            ## growing acceptor etc. lists
            new_growing_acceptor_history_list = new_growing_acceptor_history_list + current_acceptor_history
            new_growing_reaction_history_list = new_growing_reaction_history_list + current_reaction_history
            new_growing_monomer_history_list = new_growing_monomer_history_list + current_monomer_history

        growing_acceptor_history_list = new_growing_acceptor_history_list.copy()
        growing_reaction_history_list = new_growing_reaction_history_list.copy()
        growing_monomer_history_list = new_growing_monomer_history_list.copy()

        new_growing_acceptor_history_list = []
        new_growing_reaction_history_list = []
        new_growing_monomer_history_list = []

        for (a,b,c) in zip(growing_acceptor_history_list, growing_reaction_history_list, growing_monomer_history_list):
            ## get the next acceptor...
            current_acceptor_history, current_reaction_history, current_monomer_history = simulate_deprot_steps(a,b,c)

            ## growing acceptor etc. lists
            new_growing_acceptor_history_list = new_growing_acceptor_history_list + current_acceptor_history
            new_growing_reaction_history_list = new_growing_reaction_history_list + current_reaction_history
            new_growing_monomer_history_list = new_growing_monomer_history_list + current_monomer_history
        
        growing_acceptor_history_list = new_growing_acceptor_history_list.copy()
        growing_reaction_history_list = new_growing_reaction_history_list.copy()
        growing_monomer_history_list = new_growing_monomer_history_list.copy()

        new_growing_acceptor_history_list = []
        new_growing_reaction_history_list = []
        new_growing_monomer_history_list = []

    elif enzyme_of_deprot_first == 'deprot': 

        for (a,b,c) in zip(growing_acceptor_history_list, growing_reaction_history_list, growing_monomer_history_list):
            ## get the next acceptor...
            current_acceptor_history, current_reaction_history, current_monomer_history = simulate_deprot_steps(a,b,c)

            ## growing acceptor etc. lists
            new_growing_acceptor_history_list = new_growing_acceptor_history_list + current_acceptor_history
            new_growing_reaction_history_list = new_growing_reaction_history_list + current_reaction_history
            new_growing_monomer_history_list = new_growing_monomer_history_list + current_monomer_history
        
        growing_acceptor_history_list = new_growing_acceptor_history_list.copy()
        growing_reaction_history_list = new_growing_reaction_history_list.copy()
        growing_monomer_history_list = new_growing_monomer_history_list.copy()

        new_growing_acceptor_history_list = []
        new_growing_reaction_history_list = []
        new_growing_monomer_history_list = []

        #for i in range(5):
        #for i in range(3):
        for (a,b,c) in zip(growing_acceptor_history_list, growing_reaction_history_list, growing_monomer_history_list):
            ## get the next acceptor...
            current_acceptor_history, current_reaction_history, current_monomer_history = simulate_reactions_step(a,b,c, ubi_donor_list)
            ## growing acceptor etc. lists
            new_growing_acceptor_history_list = new_growing_acceptor_history_list + current_acceptor_history
            new_growing_reaction_history_list = new_growing_reaction_history_list + current_reaction_history
            new_growing_monomer_history_list = new_growing_monomer_history_list + current_monomer_history

        growing_acceptor_history_list = new_growing_acceptor_history_list.copy()
        growing_reaction_history_list = new_growing_reaction_history_list.copy()
        growing_monomer_history_list = new_growing_monomer_history_list.copy()

        new_growing_acceptor_history_list = []
        new_growing_reaction_history_list = []
        new_growing_monomer_history_list = []

    growing_acceptor_history_df = pd.DataFrame(growing_acceptor_history_list)
    growing_reaction_history_df = pd.DataFrame(growing_reaction_history_list)
    growing_monomer_history_df = pd.DataFrame(growing_monomer_history_list)

    ## if statement for acceptor, trimer, tetramer etc. etc.
    # deprotect all acceptors
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_df.map(global_deprot)
    growing_acceptor_history_max_chain_number_df = growing_acceptor_history_df.map(find_max_chain_number)
    
    ## Get the column names 
    # for index, name in enumerate(column_names): 
    column_names = list(growing_acceptor_history_df.columns)

    ## reset indexes and set types to string for whole tables
    growing_acceptor_history_df = growing_acceptor_history_df.astype(str)
    growing_reaction_history_df = growing_reaction_history_df.astype(str)
    growing_monomer_history_df = growing_monomer_history_df.astype(str)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.astype(str)
    growing_acceptor_history_max_chain_number_df = growing_acceptor_history_max_chain_number_df.astype(int)
    
    growing_acceptor_history_df = growing_acceptor_history_df.reset_index().drop('index', axis=1)
    growing_reaction_history_df = growing_reaction_history_df.reset_index().drop('index', axis=1)
    growing_monomer_history_df = growing_monomer_history_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_max_chain_number_df= growing_acceptor_history_max_chain_number_df.reset_index().drop('index', axis=1)

    print('PRE-REDUCTION 1 length: ' + str(len(growing_acceptor_history_df)))

    ### REDUCTION CODE
    ### REDUCTION FUNCTION 1: BLOCK OF CODE FOR MONOMERS THAT REACT WITH THEMSELVES
    self_reaction_table_df = pd.read_csv(".data/jsons/self_reaction_table_df.csv")
    ## here is the code that self reacts...
    self_reactions= list(self_reaction_table_df['reaction'])
    self_monomers = list(self_reaction_table_df['monomer'])
    ## USE the below... 
    for index, name in enumerate(column_names): 
        ## only for the reaction steps...
        ## kinda redundant anyway...
        ## INDEX PLAYING index = index +1 
        if enzyme_of_deprot_first == 'deprot': 
            index = index+1

        if index % 2 == 1:
            for (a, b) in zip(self_reactions, self_monomers):
                reaction_bool_df = growing_reaction_history_df[name]==str(a)
                monomer_bool_df = growing_monomer_history_df[name]==str(b)
                reaction_bool_df_copy = reaction_bool_df.copy()
                monomer_bool_df_copy = monomer_bool_df.copy()

                growing_acceptor_history_df = growing_acceptor_history_df[(monomer_bool_df & reaction_bool_df)==False]
                growing_reaction_history_df = growing_reaction_history_df[(monomer_bool_df & reaction_bool_df)==False]
                growing_monomer_history_df = growing_monomer_history_df[(monomer_bool_df & reaction_bool_df)==False]
                growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df[(monomer_bool_df & reaction_bool_df)==False]
                growing_acceptor_history_max_chain_number_df = growing_acceptor_history_max_chain_number_df[(monomer_bool_df & reaction_bool_df)==False]
                
                ## also edit the reaction and monomer df's
                reaction_bool_df = reaction_bool_df[(monomer_bool_df_copy & reaction_bool_df_copy)==False]
                monomer_bool_df = monomer_bool_df[(monomer_bool_df_copy & reaction_bool_df_copy)==False]
    
    ## reset indexes and set types to string for whole tables
    growing_acceptor_history_df = growing_acceptor_history_df.astype(str)
    growing_reaction_history_df = growing_reaction_history_df.astype(str)
    growing_monomer_history_df = growing_monomer_history_df.astype(str)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.astype(str)
    growing_acceptor_history_max_chain_number_df = growing_acceptor_history_max_chain_number_df.astype(int)
    
    growing_acceptor_history_df = growing_acceptor_history_df.reset_index().drop('index', axis=1)
    growing_reaction_history_df = growing_reaction_history_df.reset_index().drop('index', axis=1)
    growing_monomer_history_df = growing_monomer_history_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_max_chain_number_df= growing_acceptor_history_max_chain_number_df.reset_index().drop('index', axis=1)

    
    print('REDUCTION 1 length: ' + str(len(growing_acceptor_history_df)))

    ###REDUCTION FUNCTION 2: BLOCK OF CODE FOR REACTIONS THAT ADD MULTIPLE UBIS 
    max_additions = 1
    chain_number_change_df = pd.DataFrame()
    for index, name in enumerate(column_names): 
        if index == 0:
            chain_number_change_df[name]= 0
        else:
            previous_name = column_names[index-1]
            chain_number_change_df[name]= growing_acceptor_history_max_chain_number_df[name] - growing_acceptor_history_max_chain_number_df[previous_name]
        chain_number_change_df[0]= 0

    for index, name in enumerate(column_names): 
            ## only for the reaction steps...
            ## kinda redundant anyway...
            
            ## INDEX PLAYING index = index +1 
            if enzyme_of_deprot_first == 'deprot': 
                index = index+1

            if index % 2 == 1:
                growing_acceptor_history_df = growing_acceptor_history_df[chain_number_change_df[name]<=max_additions]
                growing_reaction_history_df = growing_reaction_history_df[chain_number_change_df[name]<=max_additions]
                growing_monomer_history_df = growing_monomer_history_df[chain_number_change_df[name]<=max_additions]
                growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df[chain_number_change_df[name]<=max_additions]
                growing_acceptor_history_max_chain_number_df = growing_acceptor_history_max_chain_number_df[chain_number_change_df[name]<=max_additions]
                chain_number_change_df = chain_number_change_df[chain_number_change_df[name]<=max_additions]

    ## reset indexes and set types to string for whole tables
    growing_acceptor_history_df = growing_acceptor_history_df.astype(str)
    growing_reaction_history_df = growing_reaction_history_df.astype(str)
    growing_monomer_history_df = growing_monomer_history_df.astype(str)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.astype(str)
    growing_acceptor_history_max_chain_number_df = growing_acceptor_history_max_chain_number_df.astype(int)
    
    growing_acceptor_history_df = growing_acceptor_history_df.reset_index().drop('index', axis=1)
    growing_reaction_history_df = growing_reaction_history_df.reset_index().drop('index', axis=1)
    growing_monomer_history_df = growing_monomer_history_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_max_chain_number_df= growing_acceptor_history_max_chain_number_df.reset_index().drop('index', axis=1)

    print('REDUCTION 2 length: ' + str(len(growing_acceptor_history_df)))
    ###REDUCTION FUNCTION 3: BLOCK OF CODE FOR REACTIONS THAT DONT REACT
    ### create functions that removes reactions if nothing happens during an enzyme reaction
    ### create functions that removes reactions if nothing happens with a smac reaction + this is covered by fake wash..
    reaction_didnt_happened_df = pd.DataFrame()
    for index, name in enumerate(column_names): 
        
        ## INDEX PLAYING index = index +1 
        if enzyme_of_deprot_first == 'enzyme':
            if index == 0:
                reaction_didnt_happened_df[0]= False
            ##create functions that removes reactions if nothing happens with a smac reaction        
            elif index % 2 == 0:
                previous_name = column_names[index-1]
                reaction_didnt_happened_df[name]= (growing_acceptor_history_df[name] == growing_acceptor_history_df[previous_name]) & (growing_reaction_history_df[name]=='SMAC_deprot')
            ### ### create functions that removes reactions if nothing happens....
            elif index % 2 == 1:
                previous_name = column_names[index-1]
                reaction_didnt_happened_df[name]= (growing_acceptor_history_df[name] == growing_acceptor_history_df[previous_name])

        elif enzyme_of_deprot_first == 'deprot':
            index = index +1
            if index == 1:
                reaction_didnt_happened_df[0]= False
            ##create functions that removes reactions if nothing happens with a smac reaction        
            elif index % 2 == 0:
                previous_name = column_names[index-2]
                reaction_didnt_happened_df[name]= (growing_acceptor_history_df[name] == growing_acceptor_history_df[previous_name]) & (growing_reaction_history_df[name]=='SMAC_deprot')
            ### ### create functions that removes reactions if nothing happens....
            elif index % 2 == 1:
                previous_name = column_names[index-2]
                reaction_didnt_happened_df[name]= (growing_acceptor_history_df[name] == growing_acceptor_history_df[previous_name])
            
    reaction_didnt_happened_df[0]= False
    
    for index, name in enumerate(column_names): 
        growing_acceptor_history_df = growing_acceptor_history_df[reaction_didnt_happened_df[name]==False]
        growing_reaction_history_df = growing_reaction_history_df[reaction_didnt_happened_df[name]==False]
        growing_monomer_history_df = growing_monomer_history_df[reaction_didnt_happened_df[name]==False]
        growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df[reaction_didnt_happened_df[name]==False]
        growing_acceptor_history_max_chain_number_df = growing_acceptor_history_max_chain_number_df[reaction_didnt_happened_df[name]==False]
        reaction_didnt_happened_df = reaction_didnt_happened_df[reaction_didnt_happened_df[name]==False]


    print('REDUCTION 3 length: ' + str(len(growing_acceptor_history_df)))
    ## FIGURING OUT K48 branching... 

    ## reset indexes and set types to string for whole tables
    growing_acceptor_history_df = growing_acceptor_history_df.astype(str)
    growing_reaction_history_df = growing_reaction_history_df.astype(str)
    growing_monomer_history_df = growing_monomer_history_df.astype(str)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.astype(str)
    growing_acceptor_history_max_chain_number_df = growing_acceptor_history_max_chain_number_df.astype(int)
    
    growing_acceptor_history_df = growing_acceptor_history_df.reset_index().drop('index', axis=1)
    growing_reaction_history_df = growing_reaction_history_df.reset_index().drop('index', axis=1)
    growing_monomer_history_df = growing_monomer_history_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_max_chain_number_df= growing_acceptor_history_max_chain_number_df.reset_index().drop('index', axis=1)

    for index, name in enumerate(column_names):
        if (index > 0):
            #print(index)
            ## pull out two dataframes..
            previous_name = column_names[index -1]
            current_reactant_df = growing_acceptor_history_df[previous_name][growing_reaction_history_df.loc[:,name] == 'Ube2K']
            current_product_df = growing_acceptor_history_df[name][growing_reaction_history_df.loc[:,name] == 'Ube2K']
            
            #current_reaction_df = growing_reaction_history_df[[name]][growing_reaction_history_df.loc[:,name] == 'Ube2K']
            current_reactant_product_df = pd.concat([current_reactant_df,current_product_df], axis = 1)

            current_reactant_product_df['reaction'] = list(map(find_Ube2K_or_Uce2g2, current_reactant_product_df[previous_name], current_reactant_product_df[name]))
            # map the index + values to the dataframe...
            for index, val in zip(list(current_reactant_product_df['reaction'].index), list(current_reactant_product_df['reaction'].values)):
                growing_reaction_history_df.loc[index,name] = val

    print('FIXING K48 Enyzmes 4 length: ' + str(len(growing_acceptor_history_df)))

    ## reset indexes and set types to string for whole tables
    growing_acceptor_history_df = growing_acceptor_history_df.astype(str)
    growing_reaction_history_df = growing_reaction_history_df.astype(str)
    growing_monomer_history_df = growing_monomer_history_df.astype(str)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.astype(str)
    growing_acceptor_history_max_chain_number_df = growing_acceptor_history_max_chain_number_df.astype(int)
    

    growing_acceptor_history_df = growing_acceptor_history_df.reset_index().drop('index', axis=1)
    growing_reaction_history_df = growing_reaction_history_df.reset_index().drop('index', axis=1)
    growing_monomer_history_df = growing_monomer_history_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_max_chain_number_df= growing_acceptor_history_max_chain_number_df.reset_index().drop('index', axis=1)

    ## FIGURING OUT K63 branching... 
    for index, name in enumerate(column_names):
        if (index > 0):
            #print(index)
            ## pull out two dataframes..
            previous_name = column_names[index -1]
            current_reactant_df = growing_acceptor_history_df[previous_name][growing_reaction_history_df.loc[:,name] == 'Ube13/Mms2']
            current_product_df = growing_acceptor_history_df[name][growing_reaction_history_df.loc[:,name] == 'Ube13/Mms2']
            
            #current_reaction_df = growing_reaction_history_df[[name]][growing_reaction_history_df.loc[:,name] == 'Ube2K']
            current_reactant_product_df = pd.concat([current_reactant_df,current_product_df], axis = 1)

            current_reactant_product_df['reaction'] = list(map(find_K63_branching, current_reactant_product_df[previous_name], current_reactant_product_df[name]))
            # map the index + values to the dataframe...
            for index, val in zip(list(current_reactant_product_df['reaction'].index), list(current_reactant_product_df['reaction'].values)):
                growing_reaction_history_df.loc[index,name] = val

    print('FIXING K63 Enyzmes 5 length: ' + str(len(growing_acceptor_history_df)))

    ## reset indexes and set types to string for whole tables
    growing_acceptor_history_df = growing_acceptor_history_df.astype(str)
    growing_reaction_history_df = growing_reaction_history_df.astype(str)
    growing_monomer_history_df = growing_monomer_history_df.astype(str)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.astype(str)
    growing_acceptor_history_max_chain_number_df = growing_acceptor_history_max_chain_number_df.astype(int)
    
    growing_acceptor_history_df = growing_acceptor_history_df.reset_index().drop('index', axis=1)
    growing_reaction_history_df = growing_reaction_history_df.reset_index().drop('index', axis=1)
    growing_monomer_history_df = growing_monomer_history_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_global_deprot_df = growing_acceptor_history_global_deprot_df.reset_index().drop('index', axis=1)
    growing_acceptor_history_max_chain_number_df= growing_acceptor_history_max_chain_number_df.reset_index().drop('index', axis=1)

    ## multimer size comes from 
    first_column = column_names[0]
    acceptor_size = growing_acceptor_history_max_chain_number_df[first_column].max()

    last_column = column_names[-1]
    multimer_size = growing_acceptor_history_max_chain_number_df[last_column].max()

    ### final piece of code saving everything...
    growing_acceptor_history_df.to_csv('.data/core_data/' + str(acceptor_size) + "mer__to_" + str(multimer_size) + 'mer_acceptor_history.csv', index=False)
    growing_reaction_history_df.to_csv('.data/core_data/' + str(acceptor_size) + "mer__to_" + str(multimer_size) + 'mer_final_reaction_history.csv', index=False)
    growing_monomer_history_df.to_csv('.data/core_data/' + str(acceptor_size) + "mer__to_" + str(multimer_size) + 'mer_monomer_history.csv', index=False)
    growing_acceptor_history_global_deprot_df.to_csv('.data/core_data/' + str(acceptor_size) + "mer__to_" + str(multimer_size) + 'mer_acceptor_history_global_deprot.csv', index=False)
    growing_acceptor_history_max_chain_number_df.to_csv('.data/core_data/' + str(acceptor_size) + "mer__to_" + str(multimer_size) + 'mer_acceptor_history_max_chain_number.csv', index=False)    
    
    return 