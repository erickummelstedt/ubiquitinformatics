import json 
import logging
import copy
from typing import Dict, List, Any, Union
import sys

# figure out the path issues
# home_dir = os.path.expanduser('~')
local_path = '/Users/ekummelstedt/le_code_base/ubiquitinformatics'
sys.path.insert(0, local_path)

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
Everything works; now build tests and clean up code for each of the following functions
- find_branching_site
- validate_protein_keys
- affirm_branching_sites
- validate_branching_sites
- convert_json_to_dict
- process_current_protein
- process_branch
- log_branching_details
- log_end_of_branching
- log_protein_details
- log_end_of_protein
- find_max_chain_number


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
    If any keys are missing, raises a KeyError.

    :param input_dict: Dictionary to validate
    :raises KeyError: If required keys are missing
    """
    allowed_keys = {"protein", "chain_number", "FASTA_sequence", "chain_length", "branching_sites"}

    # Check for unexpected keys
    unexpected_keys = set(input_dictionary.keys()) - allowed_keys
    if unexpected_keys:
        raise KeyError(f"Unexpected keys found: {unexpected_keys}. Allowed keys: {allowed_keys}")



### you'll want to change this so that it is not recursive but pops up in each recursive loop
def affirm_branching_sites(ubiquitin_structure):
    """
    Recursively checks if all required branching sites (M1, K6, K11, K27, K29, K33, K48, K63)
    are present at every level of the ubiquitin structure.

    If any sites are missing, it raises an AssertionError specifying the missing sites and the chain number.
    """
    REQUIRED_SITES = {"M1", "K6", "K11", "K27", "K29", "K33", "K48", "K63"}

    def _check_sites(ubiquitin_dict):
        chain_number = ubiquitin_dict["chain_number"]  # Extract chain number
        sites_in_this_ubiquitin = {site["site_name"] for site in ubiquitin_dict["branching_sites"]}

        # Check if any required sites are missing
        missing_sites = REQUIRED_SITES - sites_in_this_ubiquitin
        assert not missing_sites, (
            f"Missing sites {missing_sites} in Ubiquitin {chain_number}."
        )

    _check_sites(ubiquitin_structure)
    

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

    VALID_SITE_NAMES = {"M1", "K6", "K11", "K27", "K29", "K33", "K48", "K63"}

    for site in ubiquitin_dict.get("branching_sites", []):  # Safely get branching_sites list
        site_name = site.get("site_name", None)
        sequence_id = site.get("sequence_id", None)
        children = site.get("children", None)
        
        # Validate site_name
        if site_name not in VALID_SITE_NAMES:
            errors.append(f"Invalid site_name: {site_name}")

        # Validate sequence_id format
        if not sequence_id or "(" not in sequence_id or ")" not in sequence_id:
            errors.append(f"Invalid sequence_id format: {sequence_id}")

        # Validate children structure
        if children not in ("", "SMAC", "ABOC") and not isinstance(children, dict):
            errors.append(f"Invalid children format: {children}")

    # If errors were found, raise a single assertion error summarizing all issues
    if errors:
        raise AssertionError("\n".join(errors))
    

def convert_json_to_dict(parent_dictionary):
    """
    Converts a JSON string to a dictionary if necessary.
    If the input is already a dictionary, it remains unchanged.
    Raises a ValueError if the JSON format is invalid.

    :param input_data: JSON string or dictionary
    :return: Dictionary representation of the input data
    """

    if isinstance(parent_dictionary, str):
        try:
            # Ensure correct JSON format by replacing single quotes with double quotes
            formatted_json = parent_dictionary.replace("'", "\"")
            return json.loads(formatted_json)
        except json.JSONDecodeError as e:
            raise ValueError("Invalid JSON format: Unable to parse the string") from e
    elif isinstance(parent_dictionary, dict):
        return parent_dictionary
    else:
        raise TypeError("Input must be a dictionary or a JSON string")

    

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

## Logging functions 
def log_branching_details(branch, working_dictionary, context):
    """
    Log detailed information about the current protein and its attributes.
    """
    logging.info("===== START OF LYSINE SITE =====")
    logging.info(f"Chain Number: {working_dictionary['chain_number']}")
    logging.info(f"Lysine Site: {branch['site_name']}")
    logging.info(f"Sequence ID: {branch['sequence_id']}")
    logging.info(f"Lysine Conjugation: {branch['children']}")

def log_end_of_branching():
    """Logs the end of a lysine site branching."""
    logging.info("===== END OF LYSINE SITE =====")

def log_protein_details(working_dictionary, context):
    """
    Log detailed information about the current protein and its attributes.
    """
    logging.info(' ===== START OF PROTEIN ===== ')
    logging.info(f"Protein: {working_dictionary['protein']}")
    logging.info(f"Sequence: {working_dictionary['FASTA_sequence']}")
    logging.info(f"Chain Length List: {context['chain_length_list']}")
    logging.info(f"Chain Length: {working_dictionary['chain_length']}")
    logging.info(f"Chain Number List: {context['chain_number_list']}")
    logging.info(f"Chain Number: {working_dictionary['chain_number']}")
    logging.info(f"Branching Sites: {working_dictionary.get('branching_sites', [])}")

def log_end_of_protein(working_dictionary):
    """
    Logs the end of a protein processing step with its chain number.

    Args:
        working_dictionary (dict): The dictionary containing protein details, including 'chain_number'.
    """
    logging.info(f"===== END OF PROTEIN - UBI NUMBER: {working_dictionary['chain_number']} =====")





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
        context["multimer_string_name"] += f"his-(ubi{working_dictionary['chain_number']}"
    elif (working_dictionary["FASTA_sequence"] == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG") & (working_dictionary["chain_number"]==1):
        context["multimer_string_name"] += f"GG-(ubi{working_dictionary['chain_number']}"
    else: 
        context["multimer_string_name"] += f"(ubi{working_dictionary['chain_number']}"
    

    # Log current protein details
    log_protein_details(working_dictionary, context)

    # ensures no sites are missing
    affirm_branching_sites(working_dictionary)

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








