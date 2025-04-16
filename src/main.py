import json 
import logging
import copy
from typing import Dict, List, Any, Union
import sys

# figure out the path issues
# home_dir = os.path.expanduser('~')
# local_path = '/home/erickummelstedt/lecodebase/ubiquitinformatics/src/main.py'
local_path = '/Users/ekummelstedt/le_code_base/ubiquitinformatics'
sys.path.insert(0, local_path)

from src.utils import convert_json_to_dict
from src.logging_utils import log_protein_details, log_branching_details, log_end_of_branching, log_end_of_protein

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
From main:
- find_branching_site
- validate_protein_keys
- check_branching_sites
- check_branching_sequences
- validate_branching_sites
- check_branching_site_sequence_match
- check_children_format
- validate_branching_sites

From Utils
- match_assertion_error_contains
- all_strings_exist
- all_strings_exist_in_list
- convert_json_to_dict

Build tests for the following
Everything works; now build tests and clean up code for each of the following functions

Redo cover all: 
- iterate_through_ubiquitin
- inner_wrapper_iterate_through_ubiquitin

In utils
- process_current_protein
- process_branch

For logging_utils
- log_branching_details
- log_end_of_branching
- log_protein_details
- log_end_of_protein


NEED TO DO ThE FOLLOWING:
- iterate_through_ubiquitin (all tests on deeply nested ubiquitins)
For iterate_through_ubiquitin test the following:
- validate_protein_keys
- check_branching_sites
- check_branching_sequences
- validate_branching_sites
- check_branching_site_sequence_match
- check_children_format

Context tests left: 
- chain_number_list
- chain_length_list

Maybe add the following to iterate_through_ubiquitin
def context_template():
    return {
        "chain_number_list": [1],
        "chain_length_list": [],
        "free_lysines": [],
        "conjugated_lysines": [],
        "SMAC_lysines": [],
        "ABOC_lysines": [],
        "multimer_string_name": ""
    }

'''

def find_branching_site(sequence_id, FASTA_sequence):
    """
    Finds the branching site position in the given FASTA sequence.

    Args:
        sequence_id (str): The sequence ID containing the branching site.
        FASTA_sequence (str): The full FASTA sequence to search in.

    Returns:
        int: The 1-based index of the branching site.

    Raises:
        ValueError: If the sequence ID is not found in the FASTA sequence.
    """
    try:
        # Find the position of '(' in sequence_id to determine the amino acid index
        AA_in_sequence = sequence_id.index('(')

        # Remove parentheses from sequence ID for searching in FASTA sequence
        sequence_id_no_brackets = sequence_id.replace('(', '').replace(')', '')

        # Locate the sequence in FASTA_sequence and adjust the position
        position = FASTA_sequence.index(sequence_id_no_brackets) + AA_in_sequence + 1

        return position

    except ValueError:
        raise ValueError("substring not found")  # Ensures meaningful error handling
    

def validate_protein_keys(input_dictionary):
    """
    Validates that the given dictionary contains all required keys.
    If any keys are missing or invalid keys are present, raises a KeyError.

    :param input_dictionary: Dictionary to validate
    :raises KeyError: If required keys are missing or invalid keys are found.
    """
    allowed_keys = {"protein", "chain_number", "FASTA_sequence", "chain_length", "branching_sites"}

    # Identify missing and invalid keys
    missing_keys = allowed_keys - set(input_dictionary.keys())
    invalid_keys = set(input_dictionary.keys()) - allowed_keys

    if missing_keys or invalid_keys:
        error_messages = []
        if missing_keys:
            error_messages.append(f"Missing required keys: {missing_keys}")
        if invalid_keys:
            error_messages.append(f"Invalid keys found: {invalid_keys}")
        
        error_messages.append(f"Allowed keys: {allowed_keys}")

        raise KeyError(". ".join(error_messages))
    
def check_branching_sites(ubiquitin_dict):
    """
    Checks if all required branching sites (M1, K6, K11, K27, K29, K33, K48, K63)
    are present at that level of the ubiquitin structure.

    If any sites are missing, it raises an AssertionError specifying the missing sites and the chain number.
    """
    REQUIRED_SITES = {'M1', 'K6', 'K11', 'K27', 'K29', 'K33', 'K48', 'K63'}
    
    # sites in this ubiquitin 
    sites_in_this_ubiquitin = {site["site_name"] for site in ubiquitin_dict["branching_sites"]}
    chain_number = ubiquitin_dict["chain_number"]  # Extract chain number

    # Identify missing and invalid sites
    missing_sites = REQUIRED_SITES - sites_in_this_ubiquitin
    invalid_sites = sites_in_this_ubiquitin - REQUIRED_SITES
    
    error_messages = []
    if missing_sites or invalid_sites:
        if missing_sites:
            error_messages.append(f"Missing required sites: {missing_sites} in Ubiquitin {chain_number}")
        if invalid_sites:
            error_messages.append(f"Invalid sites found: {invalid_sites} in Ubiquitin {chain_number}")
        error_messages.append("Allowed sites: {'M1', 'K6', 'K11', 'K27', 'K29', 'K33', 'K48', 'K63'}")
    
    return error_messages

### you'll want to change this so that it is not recursive but pops up in each recursive loop
def check_branching_sequences(ubiquitin_dict):
    """
    Checks if all required branching sequence ids
    ((M)QIF, IFV(K)TLT, LTG(K)TIT, ENV(K)AKI, VKA(K)IQD, IQD(K)EGI, FAG(K)QLE, NIQ(K)EST)
    are present at that level of the ubiquitin structure.

    If any sites are missing, it raises an AssertionError specifying the missing sites and the chain number.
    """
    REQUIRED_SEQUENCE_IDS = {"(M)QIF", "IFV(K)TLT", "LTG(K)TIT", "ENV(K)AKI", "VKA(K)IQD", "IQD(K)EGI", "FAG(K)QLE", "NIQ(K)EST"}
    
    # sites in this ubiquitin 
    sequences_in_this_ubiquitin = {sequence["sequence_id"] for sequence in ubiquitin_dict["branching_sites"]}
    chain_number = ubiquitin_dict["chain_number"]  # Extract chain number
    
    # Identify missing and invalid sequences
    missing_sequences = REQUIRED_SEQUENCE_IDS - sequences_in_this_ubiquitin
    invalid_sequences = sequences_in_this_ubiquitin - REQUIRED_SEQUENCE_IDS

    # sort the sequences for testing reliability - so they always appear on the same format..

    # joining error messages 
    error_messages = []
    if missing_sequences or invalid_sequences:
        if missing_sequences:
            error_messages.append(f"Missing required sequences: {missing_sequences} in Ubiquitin {chain_number}")
        if invalid_sequences:
            error_messages.append(f"Invalid sequences found: {invalid_sequences} in Ubiquitin {chain_number}")
        error_messages.append(f"Allowed sequences: {REQUIRED_SEQUENCE_IDS}")
    
    return error_messages

def check_branching_site_sequence_match(site, errors):
    """
    Validates if the given site_name matches its expected sequence_id.

    Appends an error to `errors` if the sequence_id is incorrect.
    """
    MATCHING_SITE_SEQUENCE = {
        "M1": "(M)QIF",
        "K6": "IFV(K)TLT",
        "K11": "LTG(K)TIT",
        "K27": "ENV(K)AKI",
        "K29": "VKA(K)IQD",
        "K33": "IQD(K)EGI",
        "K48": "FAG(K)QLE",
        "K63": "NIQ(K)EST"
    }

    site_name = site.get("site_name")
    sequence_id = site.get("sequence_id")

    expected_seq = MATCHING_SITE_SEQUENCE[site_name]
    if sequence_id != expected_seq:
        errors.append(f"site_name: {site_name}, does not correspond with the sequence_id: {sequence_id}")


def check_children_format(ubiquitin_dict, site, errors):
    """
    Validates the format of the 'children' field in a branching site.

    Acceptable formats: "", "SMAC", "ABOC", or a dictionary.
    Appends an error to `errors` if invalid.
    """
  
    children = site.get("children")
    chain_number = ubiquitin_dict.get("chain_number")

    if children not in ("", "SMAC", "ABOC") and not isinstance(children, dict):
        errors.append(f"Invalid children format: {children} in Ubiquitin {chain_number}")


def validate_branching_sites(ubiquitin_dict, errors=None):
    """
    Recursively checks all branching sites in the ubiquitin dictionary 
    to ensure proper formatting.

    Collects **all** errors instead of stopping at the first failure.

    Args:
        ubiquitin_dict (dict): The ubiquitin structure to validate.
        errors (list): A list to accumulate error messages (used for recursion).

    Raises:
        AssertionError: If any errors are found in branching site formatting.
    """
    if errors is None:
        errors = []  # Initialize error list on the first call

    MATCHING_SITE_SEQUENCE = {"M1":"(M)QIF", 
                              "K6":"IFV(K)TLT", 
                              "K11":"LTG(K)TIT", 
                              "K27":"ENV(K)AKI", 
                              "K29":"VKA(K)IQD", 
                              "K33":"IQD(K)EGI", 
                              "K48":"FAG(K)QLE", 
                              "K63":"NIQ(K)EST"}
    
    site_errors = check_branching_sites(ubiquitin_dict)
    sequence_errors = check_branching_sequences(ubiquitin_dict)

    errors = errors + site_errors
    errors = errors + sequence_errors

    # If errors were found in the first two check, raise a single assertion error summarizing all issues    
    if errors:
        raise AssertionError("\n".join(errors))
    
    for site in ubiquitin_dict.get("branching_sites", []):
        check_branching_site_sequence_match(site, errors)
        check_children_format(ubiquitin_dict,site, errors)

    # If errors were found, raise a single assertion error summarizing all issues    
    if errors:
        raise AssertionError("\n".join(errors))
    

    



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
    
    # Handle protecting groups
    if branch["children"] in ["SMAC"]:
        ## add protecting group
        context["multimer_string_name"] += f"<{branch['site_name']}_SMAC>"
        context["SMAC_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name'])]]

    elif branch["children"] in ["ABOC"]:
        ## add protecting group
        context["multimer_string_name"] += f"<{branch['site_name']}_ABOC>"
        context["ABOC_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name'])]]

    # this can change later
    # Handle free lysines
    elif (branch['children'] == "") & (branch['site_name'] in ['M1','K6','K11','K27','K29','K33']): 
        print('')

    # Handle K48 & K63 lysines
    elif (branch['children'] == "") & (branch['site_name'] in ['K48', 'K63']): 
        context["free_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name'])]]

    # Handle branches that have proteins bound  
    elif isinstance(branch["children"], dict):
        context["multimer_string_name"] += f"<{branch['site_name']}_"
        branch["children"], context = inner_wrapper_iterate_through_ubiquitin(
            branch["children"], context
        )
        context["conjugated_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name'])]]
        context["multimer_string_name"] += ">"

    return branch, working_dictionary, context


def iterate_through_ubiquitin(parent_dictionary):
    
    """
    Relabel ubiquitin chain numbers, process protein branching and validate input dictionary keys.
    Generate a multimer string name based on a parent dictionary that represents proteins and ubiquitin chains.

    Args:
        parent_dictionary (dict or str): Ubiquitin structure as a dictionary or JSON string.

    Returns:
        dict: Updated polyubiquitin dictionary/JSON with relabeled values.
        context: Upadated dictionary with the information regarding the polyubiquitin. Includes;
            - chain_number_list
            - chain_length_list
            - multimer_string_name
            - max_chain_number
            - ABOC_lysines
            - SMAC_lysines
            - conjugated_lysines

    Raises:
        KeyError: If required keys are missing.
        ValueError: If the input is not a dictionary or valid JSON string.

    """
    
    # Ensure that the parent dictionary is a JSON
    parent_dictionary = convert_json_to_dict(parent_dictionary)

    # Ensure that the parent dictionary has all the valid keys 
    validate_protein_keys(parent_dictionary)
    
    # Initialize context object to track global-like variables
    context = {
        "chain_number_list": [1],
        "chain_length_list": [],
        "multimer_string_name": "",
        "max_chain_number" : "",
        "ABOC_lysines" : [],
        "SMAC_lysines": [],
        "free_lysines": [],
        "conjugated_lysines": []
    }

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

    # Ensure that the working dictionary is a JSON
    working_dictionary = convert_json_to_dict(working_dictionary)

    # Ensure that the working dictionary has all the valid keys 
    validate_protein_keys(working_dictionary)
    
    ## Process node numbers for the pre-order tree treversal of the protein 
    working_dictionary, context = process_current_protein(working_dictionary, context)

    # Append chain information to multimer string
    # Change to pdbid
    if (working_dictionary["FASTA_sequence"] == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH") & (working_dictionary["chain_number"]==1):
        context["multimer_string_name"] += f"his-GG-{working_dictionary['protein']}-{working_dictionary['chain_number']}-("
    elif (working_dictionary["FASTA_sequence"] == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG") & (working_dictionary["chain_number"]==1):
        context["multimer_string_name"] += f"GG-{working_dictionary['protein']}-{working_dictionary['chain_number']}-("
    else: 
        context["multimer_string_name"] += f"{working_dictionary['protein']}-{working_dictionary['chain_number']}-("
    
    # Log current protein details
    log_protein_details(working_dictionary, context)

    # Ensure that the branching has all the valid keys 
    validate_branching_sites(working_dictionary)

    ## get the branching sites
    working_branching_sites = working_dictionary.get("branching_sites", [])

    ## loop through the branches
    for branch in working_branching_sites:
        # log new branching deatils
        log_branching_details(branch, working_dictionary, context)
        
        ## process branch
        branch, working_dictionary, context = process_branch(branch, working_dictionary, context)

        # log end of branching details
        log_end_of_branching()

    log_end_of_protein(working_dictionary)

    # End of multimer string editing 
    context["multimer_string_name"] += ")"
    
    return working_dictionary, context

def find_max_chain_number(context):
    chain_number_list = context['chain_number_list']
    max_chain_number = chain_number_list[-1]-1
    return max_chain_number


















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

    return iterate_through_ubiquitin(adapted_dictionary)


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








