import json 
import logging
import copy
from typing import Dict, List, Any, Union
import sys

# figure out the path issues
# home_dir = os.path.expanduser('~')
# local_path = '/home/erickummelstedt/lecodebase/ubiquitinformatics/src/main.py'
local_path = '/Users/ekummelstedt/le_code_base/ubiquitinformatics/back_end'
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


















