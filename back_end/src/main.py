import json 
import logging
import copy
import sys
from pathlib import Path

# Dynamically get the backend path relative to this file
current_file = Path(__file__).resolve()
project_root = current_file.parents[2]  # Go up to project root
sys.path.insert(0, str(project_root))
local_path = project_root / 'back_end'
sys.path.insert(0, str(local_path))

from src.utils.utils import *
from src.utils.logging_utils import *

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
    working_dictionary: dict,
    context: dict
    ) -> tuple:
    
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
        context["conjugated_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name']), branch["children"]['chain_number']]]
        branch["children"], context = inner_wrapper_iterate_through_ubiquitin(
            branch["children"], context
        )
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

    # Find the maximum chain number
    context = add_max_chain_number(context)
    
    return working_dictionary, context

def add_max_chain_number(context):
    chain_number_list = context['chain_number_list']
    max_chain_number = chain_number_list[-1]-1
    context['max_chain_number'] = max_chain_number
    return context

def K_residue_ubi_addition(working_dictionary, specific_ubi_num, ubiquitination_sequence, ubi_molecule_to_add):
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
    # Ensure working_dictionary is a dict
    working_dictionary = convert_json_to_dict(working_dictionary)

    # Only convert to dict if not a protecting group
    if ubi_molecule_to_add not in ('SMAC', 'ABOC'):
        ubi_molecule_to_add = convert_json_to_dict(ubi_molecule_to_add)

    # Access immediate branching_sites
    direct_branching_sites = working_dictionary['branching_sites']

    # Loop through the lysine sites (e.g., K48, K63)
    for bra in direct_branching_sites:
        # Identify the index of the lysine within direct_branching_sites
        loop_index = list(direct_branching_sites).index(bra)

        # If lysine residue matches the ubiquitination site, proceed
        if bra['sequence_id'] == ubiquitination_sequence:
            # If the site is not conjugated (children is "", "SMAC", or "ABOC")
            if bra["children"] in ("ABOC", "SMAC", ""):
                ubi_statement = (
                    f"========== ALERT: CONJUGATION ========== this lysine {bra['site_name']} "
                    f"on ubiquitin {specific_ubi_num} is not conjugated, the ubiquitin will be added here."
                )
                logging.info(ubi_statement)
                # Ensure the C-terminus of the ubiquitin has the correct GG sequence
                if ubi_molecule_to_add['FASTA_sequence'][-5:] == 'RLRGG':
                    bra["children"] = ubi_molecule_to_add
                else:
                    raise TypeError("Ubiquitin C-terminus does not end with RLRGG.")
                # Update the direct_branching_sites with the modified branch
                direct_branching_sites[loop_index] = bra
            # Raise error if the site is already conjugated (children is a dict)
            elif isinstance(bra["children"], dict):
                raise TypeError(
                    f"The lysine specified for conjugation {bra['site_name']} on "
                    f"ubiquitin {specific_ubi_num} is already conjugated"
                )
    # Update working_dictionary with modified branching_sites
    working_dictionary['branching_sites'] = direct_branching_sites
    return working_dictionary

def ubiquitin_simulation(
        parent_dictionary: dict | str, 
        ubi_molecule_to_add: dict | str, 
        type_of_reaction: str
        ):
    """
    Simulates ubiquitin addition, deprotection, or branching reactions by recursively traversing the protein structure.

    Args:
        parent_dictionary (dict or str): The protein or ubiquitin structure in dictionary or JSON string form.
        ubi_molecule_to_add (dict or str): Ubiquitin molecule to be added, in dictionary or JSON string form.
        type_of_reaction (str): Reaction type (e.g., 'SMAC_deprot', 'GLOBAL_deprot', 'ABOC_deprot', 'K48' or 'K63').
            (type of reaction (str): is either K48 or K63, enzyme is not defined here

    Returns:
        dict: Updated dictionary representing the protein structure after the reaction.
    """

    # Validate reaction type
    valid_reactions = {'SMAC_deprot', 'GLOBAL_deprot', 'ABOC_deprot', 'FAKE_deprot', 'K48', 'K63'}
    if type_of_reaction not in valid_reactions:
        raise ValueError(f"Invalid type_of_reaction: {type_of_reaction}. Must be one of {valid_reactions}")

    # Normalize input to dicts
    parent_dictionary = convert_json_to_dict(parent_dictionary)
    if ubi_molecule_to_add not in ('SMAC', 'ABOC', ''):
        ubi_molecule_to_add = convert_json_to_dict(ubi_molecule_to_add)

    # Initialize context for chain numbering and length tracking
    context = {
        "chain_number_list": [1],
        "chain_length_list": [],
    }

    adapted_dictionary, context = inner_wrapper_ubiquitin_simulation(
        parent_dictionary, ubi_molecule_to_add, type_of_reaction, context
    )

    output_dictionary, output_context = iterate_through_ubiquitin(adapted_dictionary)

    return output_dictionary, output_context

def inner_wrapper_ubiquitin_simulation(
        input_dictionary: dict, 
        ubi_molecule_to_add: dict, 
        type_of_reaction: str,
        context
        ) -> dict:
    """
    Inner function to recursively handle ubiquitin simulation on a protein chain or structure.

    Args:
        input_dictionary (dict): Protein or ubiquitin structure for the current level of recursion.
        ubi_molecule_to_add (dict): Ubiquitin molecule to add.
        type_of_reaction (str): Type of reaction being simulated.

    Returns:
        dict: Updated protein or ubiquitin structure.
    """
    # Validate reaction type
    valid_reactions = {'SMAC_deprot', 'GLOBAL_deprot', 'ABOC_deprot', 'FAKE_deprot', 'K48', 'K63'}
    if type_of_reaction not in valid_reactions:
        raise ValueError(f"Invalid type_of_reaction: {type_of_reaction}. Must be one of {valid_reactions}")

    # Deep copy to avoid mutating input
    working_dictionary = copy.deepcopy(input_dictionary)
    working_dictionary = convert_json_to_dict(working_dictionary)

    # Set the current chain number from context
    working_dictionary['chain_number'] = context['chain_number_list'][-1]
    # Set and record chain length
    working_dictionary['chain_length'] = len(working_dictionary['FASTA_sequence'])
    context['chain_length_list'].append(working_dictionary['chain_length'])

    # Logging for debugging and traceability
    log_protein_details(working_dictionary, context)

    # Increment chain_number for future recursive calls
    context['chain_number_list'].append(context['chain_number_list'][-1] + 1)

    for bra in working_dictionary['branching_sites']:
        log_branching_details(bra, working_dictionary, context)

        # Log the type of attachment at the lysine site
        if bra['children'] in ('SMAC', 'ABOC'):
            logging.info(f"Protecting Group: {bra['children']}")
        elif bra['children'] == "":
            logging.info(f"There is no Protecting Group on: {bra['site_name']}")
        elif isinstance(bra['children'], dict):
            logging.info(f'NEXT CHAIN: {bra["children"]}')
            # Recursively process the next ubiquitin chain, updating context
            bra['children'], context = inner_wrapper_ubiquitin_simulation(
                bra['children'], ubi_molecule_to_add, type_of_reaction, context
            )

        # Actions for deprotection and conjugation logic
        bra, working_dictionary = process_ubiquitin_reaction(
            bra, working_dictionary, type_of_reaction, ubi_molecule_to_add
        )
        log_end_of_branching()

    log_end_of_protein(working_dictionary)
    return working_dictionary, context

def process_ubiquitin_reaction(
        bra: dict,
        working_dictionary: dict,
        type_of_reaction: str,
        ubi_molecule_to_add: dict
    ) -> tuple[dict, dict]:
    """
    Handle the logic for deprotection and ubiquitin conjugation at a specific lysine.

    Args:
        bra (dict): The branching site information.
        working_dictionary (dict): Current working protein or ubiquitin structure.
        type_of_reaction (str): The type of reaction being applied.
        ubi_molecule_to_add (dict): The ubiquitin molecule to add, if conjugation is triggered.

    Returns:
        tuple: Updated (bra, working_dictionary).
    """
    if type_of_reaction == 'SMAC_deprot' and bra['children'] == 'SMAC':
        bra['children'] = ''
    elif type_of_reaction == 'ABOC_deprot' and bra['children'] == 'ABOC':
        bra['children'] = ''
    elif type_of_reaction == 'GLOBAL_deprot' and bra['children'] in ('ABOC', 'SMAC'):
        bra['children'] = ''
    elif type_of_reaction == 'FAKE_deprot':
        logging.info("FAKE deprotection: No action taken.")

    # TODO: Implement other reaction types as needed
    # can add reactions for M1, K6, K11, K27, K29, K33 later
    elif (
        type_of_reaction == 'K48'
        and bra['children'] == ''
        and bra['sequence_id'] == 'FAG(K)QLE'
    ):
        working_dictionary = K_residue_ubi_addition(
            working_dictionary,
            working_dictionary['chain_number'],
            'FAG(K)QLE',
            ubi_molecule_to_add
        )
    elif (
        type_of_reaction == 'K63'
        and bra['children'] == ''
        and bra['sequence_id'] == 'NIQ(K)EST'
    ):
        working_dictionary = K_residue_ubi_addition(
            working_dictionary,
            working_dictionary['chain_number'],
            'NIQ(K)EST',
            ubi_molecule_to_add
        )
    return bra, working_dictionary




def ubiquitin_building(
    parent_dictionary: dict | str,
    ubi_molecule_to_add: dict | str,
    ubiquitin_number: int,
    lysine_residue: str
) -> dict:
    """
    Entry point for building a ubiquitin chain by attaching a molecule or protecting group.

    Args:
        parent_dictionary (dict or str): Input protein structure.
        ubi_molecule_to_add (dict or str): Ubiquitin or protecting group (SMAC/ABOC).
        ubiquitin_number (int): The chain number to which the molecule is added.
        lysine_residue (str): The specific lysine site for conjugation.

    Returns:
        dict: The updated protein dictionary with changes applied.
    """
    # Initialize context for chain numbering and length tracking
    context = {
        "chain_number_list": [1],
        "chain_length_list": [],
    }

    # Normalize input to dicts
    parent_dictionary = convert_json_to_dict(parent_dictionary)

    # Only convert to dict if not a protecting group
    if ubi_molecule_to_add not in ('SMAC', 'ABOC'):
        ubi_molecule_to_add = convert_json_to_dict(ubi_molecule_to_add)

     # Error handling fop invalid ubi_molecule_to_add
    if ubi_molecule_to_add not in ('SMAC', 'ABOC'):
        if not isinstance(ubi_molecule_to_add, dict):
            raise TypeError("ubi_molecule_to_add must be a dictionary or 'SMAC'/'ABOC' string")

    # Build structure recursively
    output_dictionary, context = inner_wrapper_ubiquitin_building(
        parent_dictionary, ubi_molecule_to_add, ubiquitin_number, lysine_residue, context
    )
    
    output_dictionary, output_context = iterate_through_ubiquitin(output_dictionary)

    return output_dictionary, output_context

def inner_wrapper_ubiquitin_building(
    input_dictionary: dict | str,
    ubi_molecule_to_add: dict | str,
    ubiquitin_number: int,
    lysine_residue: str,
    context: dict
) -> tuple:
    """
    Recursively traverse and modify the ubiquitin structure to add new branches
    or protecting groups to a given lysine site.

    Args:
        input_dictionary (dict or str): The protein or ubiquitin structure.
        ubi_molecule_to_add (dict or str): Ubiquitin or protecting group.
        ubiquitin_number (int): The specific chain number to modify.
        lysine_residue (str): Lysine site to target.
        context (dict): Tracks chain numbering and lengths.

    Returns:
        tuple: Updated working dictionary and context.
    """
    # Deep copy to avoid mutating input
    working_dictionary = copy.deepcopy(input_dictionary)
    working_dictionary = convert_json_to_dict(working_dictionary)

    # Set the current chain number from context
    working_dictionary['chain_number'] = context['chain_number_list'][-1]
    # Set and record chain length
    working_dictionary['chain_length'] = len(working_dictionary['FASTA_sequence'])
    context['chain_length_list'].append(working_dictionary['chain_length'])

    # Log protein details for debugging and traceability
    log_protein_details(working_dictionary, context)

    # Increment chain_number for future recursive calls
    context['chain_number_list'].append(context['chain_number_list'][-1] + 1)

    for bra in working_dictionary['branching_sites']:
        log_branching_details(bra, working_dictionary, context)

        # Apply modification logic to the lysine site
        bra, working_dictionary = handle_lysine_modification(
            bra, working_dictionary, ubi_molecule_to_add, ubiquitin_number, lysine_residue
        )

        # Log the state of the current branch
        if bra['children'] in ('SMAC', 'ABOC'):
            logging.info(f"Protecting Group: {bra['children']}")
        elif bra['children'] == "":
            logging.info(f"There is no Protecting Group on: {bra['site_name']}")
        elif isinstance(bra['children'], dict):
            logging.info(f"NEXT CHAIN: {bra['children']}")
            # Recursively process the next ubiquitin chain
            bra['children'], context = inner_wrapper_ubiquitin_building(
                bra['children'], ubi_molecule_to_add, ubiquitin_number, lysine_residue, context
            )
        log_end_of_branching()

    log_end_of_protein(working_dictionary)

    return working_dictionary, context

def handle_lysine_modification(
    bra: dict,
    working_dictionary: dict,
    ubi_molecule_to_add: object,
    ubiquitin_number: int,
    lysine_residue: str
) -> tuple:
    """
    Handles logic for adding a protecting group or conjugating a ubiquitin
    to a specified lysine in the target chain.

    Args:
        bra (dict): The current lysine branch site.
        working_dictionary (dict): The parent protein structure.
        ubi_molecule_to_add (str or dict): Either 'SMAC', 'ABOC', or a ubiquitin molecule.
        ubiquitin_number (int): Chain number to apply changes to.
        lysine_residue (str): Target lysine site name.

    Returns:
        tuple: (updated_branch_site, updated_working_dictionary)
    """
    if (
        ubiquitin_number == working_dictionary['chain_number'] and
        lysine_residue == bra['site_name'] and
        bra['children'] == '' and
        ubi_molecule_to_add in ('SMAC', 'ABOC')
    ):
        bra['children'] = ubi_molecule_to_add
    elif (
        ubiquitin_number == working_dictionary['chain_number'] and
        lysine_residue == bra['site_name'] and
        bra['children'] == ''
    ):
        working_dictionary = K_residue_ubi_addition(
            working_dictionary,
            working_dictionary['chain_number'],
            bra['sequence_id'],
            ubi_molecule_to_add
        )
    return bra, working_dictionary

def initialize_multimer_dicts(initial_acceptor):
    """
    Set up the initial multimer and context dictionaries.
    """
    multimer_dicts = {
        'multimers': [],
        'contexts': []
    }

    # Initialize acceptor and context
    acceptor, context = iterate_through_ubiquitin(initial_acceptor)

    # Add to the dictionary
    multimer_dicts['multimers'].append(acceptor)
    multimer_dicts['contexts'].append(context)

    return multimer_dicts


def defining_json_multimers(multimer_dicts, unprotected_ubi):
    """
    Expands a list of multimers by conjugating new ubiquitins
    at all free lysines.
    """
    multimers = multimer_dicts['multimers']
    contexts = multimer_dicts['contexts']

    new_multimer_dicts = {
        'multimers': [],
        'contexts': []
    }

    # Iterate through each multimer and its context
    for multimer, context in zip(multimers, contexts):
        free_lysines = context['free_lysines']

        for free_lysine in free_lysines:
            # Deep copy the multimer for safe modification
            working_multimer = copy.deepcopy(multimer)

            # Build new ubiquitin at the given site
            working_multimer, working_context = ubiquitin_building(
                working_multimer,
                copy.deepcopy(unprotected_ubi),
                free_lysine[0],
                free_lysine[1]
            )

            new_multimer_dicts['multimers'].append(working_multimer)
            new_multimer_dicts['contexts'].append(working_context)

    return new_multimer_dicts


def delete_duplicate_multimers(my_dict):
    """
    Removes duplicate dictionaries in each list of a given dictionary.
    """
    for key in my_dict:
        my_list = my_dict[key]
        unique = []
        seen = set()

        for d in my_list:
            marker = json.dumps(d, sort_keys=True)
            if marker not in seen:
                seen.add(marker)
                unique.append(d)

        my_dict[key] = unique

    return my_dict


