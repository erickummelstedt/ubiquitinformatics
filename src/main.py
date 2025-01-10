import json 
import logging
import copy
from typing import Dict, List, Any, Union

'''
	•	turn this into object-oriented code 
	•	correctly introduce logging 
	•	input if parent dictionary is str carry out a json loads
	•	output json loads 
	•	introduce a non-sql 
    •	fix sequence set up; M1, K6, K11, K27, K29, K33, K48, K63
'''

'''
Use descriptive names for variables, functions, and classes: Avoid abbreviations or cryptic names.
	•	Follow PEP 8 standards: The official style guide for Python code.
	•	Use 4 spaces per indentation level.
	•	Limit line length to 79 characters.
	•	Use consistent naming conventions (snake_case for functions and variables, PascalCase for classes).
'''

'''
Functions that extract information from the ubiquitin
    Pull them all together into one function --- getting_multimer_string_name 
    Only traverse the ubiquitin once
    Let it renumber 

- find_number_of_ABOC_SMAC
- find_number_of_ABOC
- find_number_of_SMAC
- relabelling_ubiquitin_numbers
- find_max_chain_number
- getting_multimer_string_name
- find_free_lysines

all these functions work off the main function

- output becomes JSON of the ubiquitin
- and the context which is a separate information file

'''

import logging
import copy
import json


def iterate_through_ubiquitin(parent_dictionary):
    """
    Relabel ubiquitin chain numbers and process protein branching.
    Generate a multimer string name based on a parent dictionary that represents proteins and ubiquitin chains.

    Args:
        parent_dictionary (Union[Dict[str, Any], str]): Dictionary representing protein data or its JSON string.
    Returns:
        Dict[str, Any]: Updated dictionary with relabeled chain numbers.
    """
    # Initialize context object to track global-like variables
    context = {
        "chain_number_list": [1],
        "chain_length_list": [],
        "multimer_string_name": ""
    }

    if isinstance(parent_dictionary, str):
        parent_dictionary = json.loads(parent_dictionary.replace("'", "\""))

    adapted_dictionary, context = inner_wrapper_iterate_through_ubiquitin(
        parent_dictionary, context
    )
    return adapted_dictionary, context

def inner_wrapper_iterate_through_ubiquitin(input_dictionary, context):
    """
    Recursively process nested dictionaries to relabel and extract ubiquitin chain numbers and process protein branching.
    Recursively process nested dictionaries to construct the multimer string name.

    """
    working_dictionary = copy.deepcopy(input_dictionary)
    
    ## Process node numbers for the pre-order tree treversal of the protein 
    working_dictionary, context = process_current_protein(working_dictionary, context)

    # Append chain information to multimer string
    # Change to pdb id
    if working_dictionary["chain_number"] == 1:
        context["multimer_string_name"] += f"his-ubi{working_dictionary['chain_number']}["
    else:
        context["multimer_string_name"] += f"ubi{working_dictionary['chain_number']}["

    # Log current protein details
    log_protein_details(working_dictionary, context)

    # Process branching sites
    working_branching_sites = working_dictionary.get("branching_sites", [])

    for branch in working_branching_sites:
        working_dictionary, context= process_branch(branch, working_dictionary, context)

    # End of multimer string editing 
    context["multimer_string_name"] += "]"
    logging.info(f"===== END OF PROTEIN - UBI NUMBER: {working_dictionary['chain_number']} =====")
    return working_dictionary, context

def process_current_protein(
    working_dictionary: Dict[str, Any],
    context: dict
    ) -> Dict[str, Any]:
    """Update chain number and chain length of current protein during recursive function."""

    # Extract the values from the context to change
    chain_number_list = context["chain_number_list"]
    chain_length_list = context["chain_length_list"]
    # Assign chain number
    working_dictionary.update({
        'chain_length': len(working_dictionary["FASTA_sequence"])
    })
    working_dictionary["chain_number"] = chain_number_list[-1]
    # Update chain length
    chain_length_list.append(working_dictionary["chain_length"])
    # Update chain number length
    chain_number_list.append(chain_number_list[-1] + 1)

    # Update the chain_number_list and the chain_length_list in the context
    context["chain_number_list"] = chain_number_list
    context["chain_length_list"] = chain_length_list
    
    return working_dictionary, context

def process_branch(branch, working_dictionary, context):
    """
    Handle the logic for each branch site in the protein structure.
    - Must deal with the following; 
    """
    chain_length_list = context["chain_length_list"]
    chain_number_list = context["chain_number_list"]
    chain_number_index = working_dictionary["chain_number"] - 1

    # Find branching site
    branching_number = find_branching_site(branch["sequence_id"], working_dictionary["FASTA_sequence"])
    ubiquitin_conjugation_site = int(sum(chain_length_list[:chain_number_index])) + branching_number

    logging.info("===== START OF LYSINE SITE =====")
    logging.info(f"Chain Number: {working_dictionary['chain_number']}")
    logging.info(f"Lysine Site: {branch['site_name']}")

    # Handle protecting groups
    if branch["children"] in ["SMAC", "ABOC"]:
        context["multimer_string_name"] += handle_protecting_group(branch)
        #context["multimer_string_name"] += f"({branch['site_name']}_"

    # Handle free lysines
    elif (branch['children'] == "") & (branch['site_name'] in ['M1','K6','K11','K27','K29','K33']): 
        logging.info(f"There is no Protecting Group on: {branch['site_name']}")

    # Handle K48 & K63 lysines
    elif (branch['children'] == "") & (branch['site_name'] in ['K48', 'K63']): 
        context["free_lysine_list"] += [working_dictionary['chain_number'], str(branch['site_name'])]
        logging.info(f"There is no Protecting Group on: {branch['site_name']}")
        #context["multimer_string_name"] += f"({branch['site_name']}_"

    # Handle branches that have proteins bound  
    elif isinstance(branch["children"], dict):
        logging.info(f"NEXT UBIQUITIN: {branch['children']}")
        context["multimer_string_name"] += f"({branch['site_name']}_"
        branch["children"], context["multimer_string_name"] = inner_wrapper_iterate_through_ubiquitin(
            branch["children"], context
        )
        context["multimer_string_name"] += ")"

    logging.info("===== END OF LYSINE SITE =====")
    return 

def handle_protecting_group(branch):
    """
    Append protecting group details to the multimer string.
    """
    protecting_group = branch["children"]
    return f"({branch['site_name']}_{protecting_group})"

def log_protein_details(working_dictionary, context):
    """
    Log detailed information about the current protein and its attributes.
    """
    logging.info(f"Protein: {working_dictionary['protein']}")
    logging.info(f"Sequence: {working_dictionary['FASTA_sequence']}")
    logging.info(f"Chain Length List: {context['chain_length_list']}")
    logging.info(f"Chain Length: {working_dictionary['chain_length']}")
    logging.info(f"Chain Number List: {context['chain_number_list']}")
    logging.info(f"Chain Number: {working_dictionary['chain_number']}")
    logging.info(f"Branching Sites: {working_dictionary.get('branching_sites', [])}")

def find_branching_site(sequence_id, FASTA_sequence):
    """
    Locate the branching site in the given FASTA sequence.
    """
    return FASTA_sequence.find(sequence_id.replace("(K)", "K"))




















'''
Functions that change the ubiquitin
- K_residue_ubi_addition
- ubiquitin_simulation
- ubiquitin_building
- global_deprot
'''


def K_residue_ubi_addition(
        working_dictionary: Union[Dict[str, Any], str], 
        specific_ubi_num: int, 
        ubiquitination_sequence: str, 
        ubi_molecule_to_add: Union[Dict[str, Any], str]) -> dict:
    """
    Adds a ubiquitin monomer to a specified lysine residue in the working dictionary.
    This function no longer adds SMAC or ABOC

    Args:
        working_dictionary (dict or str): Dictionary or JSON string representing the protein structure.
        specific_ubi_num (int): Identifier for the specific ubiquitin.
        ubiquitination_sequence (str): Sequence identifier for the lysine to be conjugated.
        ubi_molecule_to_add (dict): Ubiquitin molecule to be added.

    Returns:
        dict: Updated working dictionary with the ubiquitin molecule added.
    """
    # Parse JSON string if necessary
    if isinstance(working_dictionary, str):
        working_dictionary = json.loads(working_dictionary.replace("'", "\""))

    # Parse JSON string if necessary
    if isinstance(ubi_molecule_to_add, str):
        ubi_molecule_to_add = json.loads(ubi_molecule_to_add.replace("'", "\""))

    # Extract the branching sites
    direct_branching_sites = working_dictionary.get('branching_sites', [])
    
    for index, site in enumerate(direct_branching_sites):
        # Check if the site matches the specified ubiquitination sequence
        if site['sequence_id'] == ubiquitination_sequence:
            children = site['children']

            # If the lysine is unoccupied or occupied by a placeholder, add the new ubiquitin
            if children in {""}:
                logging.info(f"========== ALERT: CONJUGATION ==========\n"
                             f"Lysine {site['site_name']} on ubiquitin {specific_ubi_num} is not conjugated. "
                             f"The ubiquitin will be added here.")

                # Verify the ubiquitin molecule's C-terminal sequence
                if ubi_molecule_to_add['FASTA_sequence'][-5:] == 'RLRGG':
                    site['children'] = ubi_molecule_to_add
                    logging.info(f"Ubiquitin added successfully to {site['site_name']}.")
                else:
                    logging.error(f"Invalid ubiquitin molecule. Missing C-terminal 'RLRGG' sequence.")
                    return working_dictionary
            
            # Raise an error if the lysine is already conjugated to a protecting group
            elif children in ["ABOC", "SMAC"]:
                raise TypeError(f"Lysine {site['site_name']} on ubiquitin {specific_ubi_num} is protected.")

            # Raise an error if the lysine is already conjugated to another protein
            elif isinstance(children, dict):
                raise TypeError(f"Lysine {site['site_name']} on ubiquitin {specific_ubi_num} is already conjugated.")

            # Update the site in the branching sites list
            direct_branching_sites[index] = site
            break

    # Update the working dictionary with the modified branching sites
    working_dictionary['branching_sites'] = direct_branching_sites

    return working_dictionary







def ubiquitin_simulation(
        parent_dictionary: dict | str, 
        ubi_molecule_to_add: dict | str, 
        type_of_reaction: str,
        chain_number_list: List[int] = None,
        chain_length_list: List[int] = None) -> dict:
    """
    Simulates ubiquitin addition, deprotection, or branching reactions.

    Args:
        parent_dictionary (dict or str): The protein or ubiquitin structure in dictionary or JSON string form.
        ubi_molecule_to_add (dict or str): Ubiquitin molecule to be added, in dictionary or JSON string form.
        type_of_reaction (str): Reaction type (e.g., 'SMAC_deprot', 'ABOC_deprot', 'K48' or 'K63').
            (type of reaction (str): is either K48 or K63, enzyme is not defined here
        chain_number_list (List[int], optional): Tracks chain numbers. Defaults to [1].
        chain_length_list (List[int], optional): Tracks chain lengths. Defaults to [].

    Returns:
        dict: Updated dictionary representing the protein structure after the reaction.
    """
    # Assign values to chain_number_list and chain_length_list
    if chain_number_list is None:
        chain_number_list = [1]
    if chain_length_list is None:
        chain_length_list = []

    # Convert JSON strings to dictionaries if necessary
    if isinstance(parent_dictionary, str):
        parent_dictionary = json.loads(parent_dictionary.replace("'", "\""))

    if isinstance(ubi_molecule_to_add, str) and ubi_molecule_to_add != '':
        ubi_molecule_to_add = json.loads(ubi_molecule_to_add.replace("'", "\""))

    # Perform the inner simulation and re-label ubiquitin numbers
    adapted_dictionary = inner_wrapper_ubiquitin_simulation(
        parent_dictionary, ubi_molecule_to_add, type_of_reaction, chain_number_list, chain_length_list
        )

    return relabelling_ubiquitin_numbers(adapted_dictionary)


def inner_wrapper_ubiquitin_simulation(
        input_dictionary: dict, 
        ubi_molecule_to_add: dict, 
        type_of_reaction: str,
        chain_number_list: List[int],
        chain_length_list: List[int]
        ) -> dict:
    """
    Inner function to recursively handle ubiquitin simulation on a protein chain or structure.

    Args:
        input_dictionary (dict): Protein or ubiquitin structure for the current level of recursion.
        ubi_molecule_to_add (dict): Ubiquitin molecule to add.
        type_of_reaction (str): Type of reaction being simulated.
        chain_number_list (List[int], optional): Tracks chain numbers. Defaults to [1].
        chain_length_list (List[int], optional): Tracks chain lengths. Defaults to [].

    Returns:
        dict: Updated protein or ubiquitin structure.
    """
    # Deep copy the dictionary to avoid mutating the original
    working_dictionary = copy.deepcopy(input_dictionary)

    # Update chain number and chain length
    working_dictionary['chain_number'] = chain_number_list[-1]
    working_dictionary['chain_length'] = len(working_dictionary['FASTA_sequence'])
    chain_length_list.append(working_dictionary['chain_length'])

    # Logging key information
    log_protein_details(working_dictionary)

    # Update chain number for subsequent branches
    chain_number_list.append(chain_number_list[-1] + 1)

    # Loop through branching sites for reactions
    working_branching_sites = working_dictionary['branching_sites']
    for bra in working_branching_sites:
        logging.info("===== START OF LYSINE SITE =====")
        logging.info(f"Chain Number: {working_dictionary['chain_number']}")
        logging.info(f"Lysine Site: {bra['site_name']}")

        # Determine the current state of the lysine site
        if bra['children'] in ['SMAC', 'ABOC']:
            logging.info(f"Protecting Group: {bra['children']}")
        elif bra['children'] == "":
            logging.info(f"Lysine {bra['site_name']} is unoccupied.")
        elif isinstance(bra['children'], dict):
            # Recursively process nested ubiquitin chains
            logging.info(f"Next Chain: {bra['children']}")
            bra['children'] = inner_wrapper_ubiquitin_simulation(
                bra['children'], ubi_molecule_to_add, type_of_reaction, chain_number_list, chain_length_list
                )

        # Perform actions based on the type of reaction
        if type_of_reaction == 'SMAC_deprot' and bra['children'] == 'SMAC':
            bra['children'] = ''
        elif type_of_reaction == 'ABOC_deprot' and bra['children'] == 'ABOC':
            bra['children'] = ''
        elif type_of_reaction == 'GLOBAL_deprot' and bra['children'] in ['ABOC', 'SMAC']:
            bra['children'] = ''
        elif type_of_reaction == 'K48' and bra['children'] == '' and bra['sequence_id'] == 'FAG(K)QLE':
            working_dictionary = K_residue_ubi_addition(working_dictionary, working_dictionary['chain_number'], 'FAG(K)QLE', ubi_molecule_to_add)
        elif type_of_reaction == 'K63' and bra['children'] == '' and bra['sequence_id'] == 'NIQ(K)EST':
            working_dictionary = K_residue_ubi_addition(working_dictionary, working_dictionary['chain_number'], 'NIQ(K)EST', ubi_molecule_to_add)

        # Log the end of lysine site processing
        logging.info("===== END OF LYSINE SITE =====")

    # Log end of chain processing
    logging.info(f"===== END OF PROTEIN - CHAIN NUMBER: {working_dictionary['chain_number']} =====")
    return working_dictionary








