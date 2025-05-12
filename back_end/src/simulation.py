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
    inject_branching_sites

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



def find_Ube2K_or_Uce2g2(reactant_dictionary, product_dictionary):
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
    adapted_dictionary = inner_wrapper_find_Ube2K_or_Uce2g2(reactant_dictionary)

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
    adapted_dictionary = inner_wrapper_find_Ube2K_or_Uce2g2(product_dictionary)

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

    #if not(specific_ubi_num in chain_number_list[:-1]): $
    #    raise TypeError('the ubiquitin specific is outside the range of the ubiquitin multimer')

### go into the protein, then the branching_sites... then loop the problem... 
### this function loops through the ubiquitin labeling the proteins...
def inner_wrapper_find_Ube2K_or_Uce2g2(input_dictionary):
    ### super important because the dictionary is a mutable object
    working_dictionary = copy.deepcopy(input_dictionary)

    ## set the chain_number from the chain_number_list
    working_dictionary['chain_number'] = int(chain_number_list[-1])
    
    ### add the current chain length to the chain length list this is useful for later for numbering of amino acids
    working_dictionary['chain_length'] = len(working_dictionary['FASTA_sequence'])
    chain_length_list.append(working_dictionary['chain_length'])

    ## printing information 
    logging.info("Protein: " + str(working_dictionary['protein']))
    logging.info("Sequence: " + str(working_dictionary['FASTA_sequence']))
    logging.info("Chain Length List: " + str(chain_number_list))
    logging.info("Chain Length: " + str(working_dictionary['chain_length']))
    logging.info("Chain Number List: " + str(chain_number_list))
    logging.info("Chain Number: " + str(working_dictionary['chain_number']))
    logging.info("Branching Sites: " + str(working_dictionary['branching_sites']))

    ### provide the branching_sites of the protein to the working_branching_sites variable 
    working_branching_sites = working_dictionary['branching_sites']
        
    ### now that we have associated a value to the current chain, we can add to the chain list
    chain_number_list.append((chain_number_list[-1]+1))
        
    ### working_branching_sites = loop_through_immediate_branching_sites(working_branching_sites, x_ubi_num)
    ### working_dictionary['branching_sites'] = working_branching_sites
    for bra in working_branching_sites:
        
        logging.info(' ===== START OF LYSINE SITE =====  ')
        logging.info("Chain Number: " + str(working_dictionary['chain_number']))
        ### which lysine site.... 
        logging.info("Branching Site: " + str(bra['site_name'])) 

        ## some print work
        if (bra['children'] =='SMAC' or bra['children'] =='ABOC'):
            logging.info("Protecting Group: " + str(bra['children']))
        
        ### if the site is a protein.... 
        elif (bra['children'] == "") & (bra['site_name'] in ['K11','K6','K27','K29','K33','M1']): 
            logging.info("There is no Protecting Group on: " + str(bra['site_name']))

        elif (bra['children'] == "") & (bra['site_name'] in ['K48', 'K63']): 
            free_lysine_list.append([working_dictionary['chain_number'], str(bra['site_name'])])
            logging.info("There is no Protecting Group on: " + str(bra['site_name']))    
        
        elif isinstance(bra['children'], dict): 
            ### recursive function that calls itself... 
            logging.info('NEXT CHAIN: ' + str(bra['children']))
            bound_lysine_list.append([working_dictionary['chain_number'], str(bra['site_name'])])
            ### the difference.... 
            bra['children'] = inner_wrapper_find_Ube2K_or_Uce2g2(bra['children']) 

        ### adding branching_sites... here can use the chain_number and lysine site...
        ### ADD ERROR for ubiquitin not found... 
        #if ubiquitin_number == working_dictionary['chain_number'] and lysine_number == None and protecting_group == None:
        #    working_dictionary['branching_sites'] = [{"site": "K48", "children": ""},
        #                                        {"site": "K11", "children": ""},
        #                                        {"site": "K63", "children": ""},
        #                                        {"site": "K6", "children": ""},
        #                                        {"site": "K27", "children": ""},
        #                                        {"site": "K29", "children": ""},
        #                                        {"site": "K33", "children": ""},
        #                                        {"site": "M1", "children": ""}]

        
       ### ADD ERROR an elif or else... which should be an error... for protecting group nor ubiquitin found... 
        logging.info(' ===== END OF LYSINE SITE =====  ')
    logging.info(' ===== END OF PROTEIN - CHAIN NUMBER: ' + str(working_dictionary['chain_number']) + ' =====  ')
    ## finish with diving into ubiquitin to relabel or the ubiquitins...
    return working_dictionary



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

### go into the protein, then the branching_sites... then loop the problem... 
### this function loops through the ubiquitin labeling the proteins...
def inner_wrapper_find_K63_branching(input_dictionary):
    ### super important because the dictionary is a mutable object
    working_dictionary = copy.deepcopy(input_dictionary)

    ## set the chain_number from the chain_number_list
    working_dictionary['chain_number'] = int(chain_number_list[-1])
    
    ### add the current chain length to the chain length list this is useful for later for numbering of amino acids
    working_dictionary['chain_length'] = len(working_dictionary['FASTA_sequence'])
    chain_length_list.append(working_dictionary['chain_length'])

    ## printing information 
    logging.info("Protein: " + str(working_dictionary['protein']))
    logging.info("Sequence: " + str(working_dictionary['FASTA_sequence']))
    logging.info("Chain Length List: " + str(chain_number_list))
    logging.info("Chain Length: " + str(working_dictionary['chain_length']))
    logging.info("Chain Number List: " + str(chain_number_list))
    logging.info("Chain Number: " + str(working_dictionary['chain_number']))
    logging.info("Branching Sites: " + str(working_dictionary['branching_sites']))

    ### provide the branching_sites of the protein to the working_branching_sites variable 
    working_branching_sites = working_dictionary['branching_sites']
        
    ### now that we have associated a value to the current chain, we can add to the chain list
    chain_number_list.append((chain_number_list[-1]+1))
        
    ### working_branching_sites = loop_through_immediate_branching_sites(working_branching_sites, x_ubi_num)
    ### working_dictionary['branching_sites'] = working_branching_sites
    for bra in working_branching_sites:
        
        logging.info(' ===== START OF LYSINE SITE =====  ')
        logging.info("Chain Number: " + str(working_dictionary['chain_number']))
        ### which lysine site.... 
        logging.info("Branching Site: " + str(bra['site_name'])) 

        ## some print work
        if (bra['children'] =='SMAC' or bra['children'] =='ABOC'):
            logging.info("Protecting Group: " + str(bra['children']))
        
        ### if the site is a protein.... 
        elif (bra['children'] == "") & (bra['site_name'] in ['K11','K6','K27','K29','K33','M1']): 
            logging.info("There is no Protecting Group on: " + str(bra['site_name']))

        elif (bra['children'] == "") & (bra['site_name'] in ['K48', 'K63']): 
            free_lysine_list.append([working_dictionary['chain_number'], str(bra['site_name'])])
            logging.info("There is no Protecting Group on: " + str(bra['site_name']))    
        
        elif isinstance(bra['children'], dict): 
            ### recursive function that calls itself... 
            logging.info('NEXT CHAIN: ' + str(bra['children']))
            bound_lysine_list.append([working_dictionary['chain_number'], str(bra['site_name'])])
            ### the difference.... 
            bra['children'] = inner_wrapper_find_K63_branching(bra['children']) 

        ### adding branching_sites... here can use the chain_number and lysine site...
        ### ADD ERROR for ubiquitin not found... 
        #if ubiquitin_number == working_dictionary['chain_number'] and lysine_number == None and protecting_group == None:
        #    working_dictionary['branching_sites'] = [{"site": "K48", "children": ""},
        #                                        {"site": "K11", "children": ""},
        #                                        {"site": "K63", "children": ""},
        #                                        {"site": "K6", "children": ""},
        #                                        {"site": "K27", "children": ""},
        #                                        {"site": "K29", "children": ""},
        #                                        {"site": "K33", "children": ""},
        #                                        {"site": "M1", "children": ""}]

        
       ### ADD ERROR an elif or else... which should be an error... for protecting group nor ubiquitin found... 
        logging.info(' ===== END OF LYSINE SITE =====  ')
    logging.info(' ===== END OF PROTEIN - CHAIN NUMBER: ' + str(working_dictionary['chain_number']) + ' =====  ')
    ## finish with diving into ubiquitin to relabel or the ubiquitins...
    return working_dictionary



