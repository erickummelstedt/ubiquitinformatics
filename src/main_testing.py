import json 
import logging
import copy
from typing import Dict, List, Any, Union
import sys

'''
Functions that extract information from the ubiquitin
    Pull them all together into one function --- getting_multimer_string_name 
    Only traverse the ubiquitin once
    Let it renumber 

- DONE
- find_branching_site
- validate_branching_sites
- relabelling_ubiquitin_numbers
- inner_wrapper_relabelling_ubiquitin_numbers
- getting_multimer_string_name
- find_number_of_ABOC_SMAC
- find_number_of_ABOC
- find_number_of_SMAC
- find_max_chain_number
- validate all branching sites



- find_free_lysines

all these functions work off the main function

- output becomes JSON of the ubiquitin
- and the context which is a separate information file

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
    
import json

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

    VALID_SITE_NAMES = {"K6", "K11", "K27", "K29", "K33", "K48", "K63", "M1"}

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

        # Recursively validate nested children
        if isinstance(children, dict):
            validate_branching_sites(children, errors)

    # If errors were found, raise a single assertion error summarizing all issues
    if errors:
        raise AssertionError("\n".join(errors))


def relabelling_ubiquitin_numbers(parent_dictionary):
    """
    Function to relabel ubiquitin numbers and validate input dictionary keys.

    Args:
        parent_dictionary (dict or str): Ubiquitin structure as a dictionary or JSON string.

    Returns:
        dict: Updated ubiquitin dictionary with relabeled values.

    Raises:
        KeyError: If required keys are missing.
        ValueError: If the input is not a dictionary or valid JSON string.
    """
    # Define the required keys
    required_keys = {"protein", "chain_number", "FASTA_sequence", "chain_length", "branching_sites"}

    global chain_number_list
    global chain_length_list
    chain_number_list = [1]
    chain_length_list = []

    # Convert JSON string to dictionary if necessary
    if isinstance(parent_dictionary, str):
        try:
            x = parent_dictionary.replace("'", "\"")  # Ensure correct JSON format
            parent_dictionary = json.loads(x)
        except json.JSONDecodeError as e:
            raise ValueError("Invalid JSON format") from e

    # Ensure input is a dictionary
    if not isinstance(parent_dictionary, dict):
        raise ValueError("Input must be a dictionary or a valid JSON string")

    # Check for missing required keys
    missing_keys = required_keys - set(parent_dictionary.keys())
    if missing_keys:
        raise KeyError(f"Missing required keys: {missing_keys}")

    # Retrieve values from the validated dictionary
    next_working_protein = parent_dictionary['protein']
    next_working_FASTA_sequence = parent_dictionary['FASTA_sequence']

    # Process the dictionary using the inner function
    adapted_dictionary = inner_wrapper_relabelling_ubiquitin_numbers(parent_dictionary)

    return adapted_dictionary


### go into the protein, then the branching_sites... then loop the problem... 
### this function loops through the ubiquitin labeling the proteins...
def inner_wrapper_relabelling_ubiquitin_numbers(input_dictionary):  
    
    # Allowed keys
    allowed_keys = {"protein", "chain_number", "FASTA_sequence", "chain_length", "branching_sites"}

    # Ensure input is a dictionary
    if not isinstance(input_dictionary, dict):
        raise ValueError("Input must be a dictionary or a valid JSON string")

    # Check for unexpected keys
    unexpected_keys = set(input_dictionary.keys()) - allowed_keys
    if unexpected_keys:
        raise KeyError(f"Unexpected keys found: {unexpected_keys}. Allowed keys: {allowed_keys}")

    # Check for missing required keys
    missing_keys = allowed_keys - set(input_dictionary.keys())
    if missing_keys:
        raise KeyError(f"Missing required keys: {missing_keys}")
    
    ### super important because the dictionary is a mutable object
    working_dictionary = copy.deepcopy(input_dictionary)
    #working_ubiAllAtoms = input_ubiAllAtoms
    #working_ubiAllBonds = input_ubiAllBonds
     
    logging.info(' ===== START OF PROTEIN ===== ')

    ## set the chain_number from the chain_number_list
    working_dictionary['chain_number'] = chain_number_list[-1]
    
    ### add the current chain length to the chain length list this is useful for later for numbering of amino acids
    
    working_dictionary['chain_length'] = len(working_dictionary['FASTA_sequence'])
    chain_length_list.append(working_dictionary['chain_length'])

    ## printing information 
    logging.info("Protein: " + str(working_dictionary['protein']))
    logging.info("Sequence: " + str(working_dictionary['FASTA_sequence']))
    logging.info("Chain Length List: " + str(chain_length_list))
    logging.info("Chain Length: " + str(working_dictionary['chain_length']))
    logging.info("Chain Number List: " + str(chain_number_list))
    logging.info("Chain Number: " + str(working_dictionary['chain_number']))
    logging.info("Branching Sites: " + str(working_dictionary['branching_sites']))

    ### provide the branching_sites of the protein to the working_branching_sites variable 
    working_branching_sites = working_dictionary['branching_sites']
    working_FASTA_sequence = working_dictionary['FASTA_sequence']
        
    ### now that we have associated a value to the current chain, we can add to the chain list
    chain_number_list.append((chain_number_list[-1]+1))

    ### adding branching_sites... here can use the chain_number and lysine site...
    ### working_branching_sites = loop_through_immediate_branching_sites(working_branching_sites, x_ubi_num)
    ### working_dictionary['branching_sites'] = working_branching_sites
    for bra in working_branching_sites:
        
        ## the following two lines finds branching site.... of the particular protein of interest
        working_sequence_id = bra['sequence_id']
        branching_number = find_branching_site(working_sequence_id, working_FASTA_sequence)

        ## add in the addition of the total chain length
        ## included is the -1 because with sum the index starts at 0
        chain_number_index = working_dictionary['chain_number']-1
        ## 
        ubiquitinConjugationSite = int(sum(chain_length_list[:chain_number_index])) + branching_number
        conjugatingAtomsStart = int(sum(chain_length_list[:chain_number_index])) + 1
        conjugatingAtomsEnd = int(sum(chain_length_list[:(chain_number_index+1)]))
        lastBackBoneSubstID = int(sum(chain_length_list)) 

        logging.info(' ===== START OF LYSINE SITE =====  ')
        logging.info("Chain Number: " + str(working_dictionary['chain_number']))
        ### which lysine site.... 
        logging.info("Lysine Site: " + str(bra['site_name']))

        ### what is attached to the lysine
        ### three options; 
        ## there is a protecting group 
        ## the lysine is free
        ## another protein is attached...
        if (bra['children'] =='SMAC' or bra['children'] =='ABOC'):
            logging.info("Protecting Group: " + str(bra['children']))
            if bra['children'] =='SMAC':
                protecting_group = 'SMAC'
            elif bra['children'] =='ABOC':
                protecting_group = 'ABOC'
                
            ### download new ubiquitin molecule
            #if protecting_group == 'SMAC':
            #    new_PG_df = pd.read_csv(local_path + ".data/mol2_files/protecting_groups/SMAC_Cl.mol2", names=['#'])
            #elif protecting_group == 'ABOC':
            #    new_PG_df = pd.read_csv(local_path + ".data/mol2_files/protecting_groups/ABOC_Cl.mol2", names=['#'])



            #new_AllAtoms, new_AllBonds = get_BONDS_ATOMS(new_PG_df)
            ### rotate new ubiquitin molecule
            #newAtoms, conjugatingAtoms, newBonds, conjugatingUbiBonds = rotation_transformation_of_new_PG(working_ubiAllAtoms, working_ubiAllBonds, new_AllAtoms, new_AllBonds, ubiquitinConjugationSite, conjugatingAtomsStart, conjugatingAtomsEnd, protecting_group)
            ### renumber the bonding for the ubiquitin molecule
            ### 

            #working_ubiAllAtoms, working_ubiAllBonds = numbering_for_bonds_and_atoms_PG(newAtoms, conjugatingAtoms, newBonds, conjugatingUbiBonds, ubiquitinConjugationSite, protecting_group)
            #logging.info(' ===== END OF PROTECTING GROUP SITE =====  ')



        ### if the site is a protein.... 
        elif bra['children'] == "": 
            logging.info("There is no Protecting Group on: " + str(bra['site_name']))
        
        elif isinstance(bra['children'], dict): 
            ### recursive function that calls itself... 
            logging.info('NEXT UBIQUITIN: ' + str(bra['children']))
            ## add the next ubiquitin here
            ## all the changes occur here
            
            ## find conjugation site
            
            
            ## NEXT WORKING PROTEIN FUNCTION ===================  NEED TO DO
            ## dive into the next dictionary.... 
            next_working_protein = bra['children']['protein']
            next_working_FASTA_sequence = bra['children']['FASTA_sequence']
            ### check the sequence...
            ### from the sequence pull the correct ubiquitin mol2... 
            ### new_ubi_df becomes the correct ubiquitin mol2...
            #if next_working_FASTA_sequence == 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG':
            #    new_ubi_df = pd.read_csv(local_path + ".data/mol2_files/monomers/single_ubiquitin.txt", names=['#'])
            #elif next_working_FASTA_sequence == 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH':
            #    new_ubi_df = pd.read_csv(local_path + ".data/mol2_files/monomers/histag_ubiquitin.txt", names=['#'])



            ### download the correct ubiquiitin molecule
            #new_ubiAllAtoms, new_ubiAllBonds = get_BONDS_ATOMS(new_ubi_df)
            ### rotate new ubiquitin molecule
            #newUbiAtoms, conjugatingAtoms, newUbiBonds, conjugatingUbiBonds = rotation_transformation_of_new_ubiquitin(working_ubiAllAtoms, working_ubiAllBonds, new_ubiAllAtoms, new_ubiAllBonds, ubiquitinConjugationSite, conjugatingAtomsStart, conjugatingAtomsEnd)
            ### renumber the bonding for the ubiquitin molecule
            #concatUbiAtoms, concatUbiBonds = numbering_for_bonds_and_atoms(newUbiAtoms, conjugatingAtoms, newUbiBonds, conjugatingUbiBonds, ubiquitinConjugationSite, lastBackBoneSubstID)
            #


            ## ====================
            ## new ubiquitin addition in JSON, and in 
            #bra['children'], working_ubiAllAtoms, working_ubiAllBonds = inner_wrapper_relabelling_ubiquitin_numbers_and_JSON_to_MOL2(bra['children'], concatUbiAtoms, concatUbiBonds) 
            bra['children'] = inner_wrapper_relabelling_ubiquitin_numbers(bra['children']) 
            
        logging.info(' ===== END OF LYSINE SITE =====  ')

        logging.info(' ===== END OF LYSINE SITE =====  ')

    logging.info(' ===== END OF PROTEIN - UBI NUMBER: ' + str(working_dictionary['chain_number']) + ' =====  ')

    return working_dictionary

def getting_multimer_string_name(parent_dictionary):
    global chain_number_list
    global chain_length_list
    global multimer_string_name
    chain_number_list = [1]
    chain_length_list = []
    multimer_string_name = ''    

    ## here download the original ubiquitin...
    ## here check which ubiquitin has is your starting ubiquitin and download that one... 
    
    if isinstance(parent_dictionary, str):
        x = parent_dictionary.replace("'", "\"")
        parent_dictionary = json.loads(x)

    ## NEXT WORKING PROTEIN FUNCTION ===================  NEED TO DO
    ## dive into the next dictionary.... 
    next_working_protein = parent_dictionary['protein']
    next_working_FASTA_sequence = parent_dictionary['FASTA_sequence']
    ## check the sequence...
    ## from the sequence pull the correct ubiquitin mol2... 
    ## new_ubi_df becomes the correct ubiquitin mol2...
    #if next_working_FASTA_sequence == 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG':
    #    start_ubi_df = pd.read_csv(local_path + ".data/mol2_files/monomers/single_ubiquitin.txt", names=['#'])
    #elif next_working_FASTA_sequence == 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH':
    #    start_ubi_df = pd.read_csv(local_path + ".data/mol2_files/monomers/histag_ubiquitin.txt", names=['#'])
    #start_ubiAllAtoms, start_ubiAllBonds = get_BONDS_ATOMS(start_ubi_df)

    #adapted_dictionary, final_ubiAllAtoms, final_ubiAllBonds = inner_wrapper_relabelling_ubiquitin_numbers_and_JSON_to_MOL2(parent_dictionary,start_ubiAllAtoms, start_ubiAllBonds)
    #write_mol2(final_ubiAllAtoms,final_ubiAllBonds)
    
    adapted_dictionary, multimer_string_name = inner_wrapper_getting_multimer_string_name(parent_dictionary, multimer_string_name)
    return multimer_string_name

### go into the protein, then the branching_sites... then loop the problem... 
### this function loops through the ubiquitin labeling the proteins...
def inner_wrapper_getting_multimer_string_name(input_dictionary, multimer_string_name):  
    ### super important because the dictionary is a mutable object
    working_dictionary = copy.deepcopy(input_dictionary)
    #working_ubiAllAtoms = input_ubiAllAtoms
    #working_ubiAllBonds = input_ubiAllBonds
     
    logging.info(' ===== START OF PROTEIN ===== ')

    ## set the chain_number from the chain_number_list
    working_dictionary['chain_number'] = chain_number_list[-1]
    if working_dictionary['chain_number'] == 1: 
        multimer_string_name = multimer_string_name + 'his-ubi' + str(working_dictionary['chain_number']) +'['
    else: 
        multimer_string_name = multimer_string_name + 'ubi' + str(working_dictionary['chain_number']) +'['

    ### add the current chain length to the chain length list this is useful for later for numbering of amino acids
    
    working_dictionary['chain_length'] = len(working_dictionary['FASTA_sequence'])
    chain_length_list.append(working_dictionary['chain_length'])

    ## printing information 
    logging.info("Protein: " + str(working_dictionary['protein']))
    logging.info("Sequence: " + str(working_dictionary['FASTA_sequence']))
    logging.info("Chain Length List: " + str(chain_length_list))
    logging.info("Chain Length: " + str(working_dictionary['chain_length']))
    logging.info("Chain Number List: " + str(chain_number_list))
    logging.info("Chain Number: " + str(working_dictionary['chain_number']))
    logging.info("Branching Sites: " + str(working_dictionary['branching_sites']))

    ### provide the branching_sites of the protein to the working_branching_sites variable 
    working_branching_sites = working_dictionary['branching_sites']
    working_FASTA_sequence = working_dictionary['FASTA_sequence']
        
    ### now that we have associated a value to the current chain, we can add to the chain list
    chain_number_list.append((chain_number_list[-1]+1))

    ### adding branching_sites... here can use the chain_number and lysine site...
    ### working_branching_sites = loop_through_immediate_branching_sites(working_branching_sites, x_ubi_num)
    ### working_dictionary['branching_sites'] = working_branching_sites
    first_dash = True

    for bra in working_branching_sites:

        ## the following two lines finds branching site.... of the particular protein of interest
        working_sequence_id = bra['sequence_id']
        branching_number = find_branching_site(working_sequence_id, working_FASTA_sequence)

        ## add in the addition of the total chain length
        ## included is the -1 because with sum the index starts at 0
        chain_number_index = working_dictionary['chain_number']-1
        ## 
        ubiquitinConjugationSite = int(sum(chain_length_list[:chain_number_index])) + branching_number
        conjugatingAtomsStart = int(sum(chain_length_list[:chain_number_index])) + 1
        conjugatingAtomsEnd = int(sum(chain_length_list[:(chain_number_index+1)]))
        lastBackBoneSubstID = int(sum(chain_length_list)) 

        logging.info(' ===== START OF LYSINE SITE =====  ')
        logging.info("Chain Number: " + str(working_dictionary['chain_number']))
        ### which lysine site.... 
        logging.info("Lysine Site: " + str(bra['site_name']))

        ### what is attached to the lysine
        ### three options; 
        ## there is a protecting group 

        ## the lysine is free
        ## another protein is attached...
        if (bra['children'] =='SMAC' or bra['children'] =='ABOC'):
            logging.info("Protecting Group: " + str(bra['children']))
            if bra['children'] =='SMAC':
                protecting_group = 'SMAC'
                if first_dash:
                    first_dash = False
                    multimer_string_name = multimer_string_name + str(bra['site_name']) + '_'
                else: 
                    multimer_string_name = multimer_string_name + '-' + str(bra['site_name']) + '_'
                multimer_string_name = multimer_string_name + 'SMAC'

            elif bra['children'] =='ABOC':
                protecting_group = 'ABOC'
                if first_dash:
                    first_dash = False
                    multimer_string_name = multimer_string_name + str(bra['site_name']) + '_'
                else: 
                    multimer_string_name = multimer_string_name + '-' + str(bra['site_name']) + '_'
                multimer_string_name = multimer_string_name + 'ABOC'
                
            ### download new ubiquitin molecule
            #if protecting_group == 'SMAC':
            #    new_PG_df = pd.read_csv(local_path + ".data/mol2_files/protecting_groups/SMAC_Cl.mol2", names=['#'])
            #elif protecting_group == 'ABOC':
            #    new_PG_df = pd.read_csv(local_path + ".data/mol2_files/protecting_groups/ABOC_Cl.mol2", names=['#'])



            #new_AllAtoms, new_AllBonds = get_BONDS_ATOMS(new_PG_df)
            ### rotate new ubiquitin molecule
            #newAtoms, conjugatingAtoms, newBonds, conjugatingUbiBonds = rotation_transformation_of_new_PG(working_ubiAllAtoms, working_ubiAllBonds, new_AllAtoms, new_AllBonds, ubiquitinConjugationSite, conjugatingAtomsStart, conjugatingAtomsEnd, protecting_group)
            ### renumber the bonding for the ubiquitin molecule
            ### 

            #working_ubiAllAtoms, working_ubiAllBonds = numbering_for_bonds_and_atoms_PG(newAtoms, conjugatingAtoms, newBonds, conjugatingUbiBonds, ubiquitinConjugationSite, protecting_group)
            #logging.info(' ===== END OF PROTECTING GROUP SITE =====  ')



        ### if the site is a protein.... 
        elif bra['children'] == "": 
            logging.info("There is no Protecting Group on: " + str(bra['site_name']))
        
        elif isinstance(bra['children'], dict): 
            ### recursive function that calls itself... 
            logging.info('NEXT UBIQUITIN: ' + str(bra['children']))
            ## add the next ubiquitin here
            ## all the changes occur here
            
            ## find conjugation site
            
            
            ## NEXT WORKING PROTEIN FUNCTION ===================  NEED TO DO
            ## dive into the next dictionary.... 
            next_working_protein = bra['children']['protein']
            next_working_FASTA_sequence = bra['children']['FASTA_sequence']
            ### check the sequence...
            ### from the sequence pull the correct ubiquitin mol2... 
            ### new_ubi_df becomes the correct ubiquitin mol2...
            #if next_working_FASTA_sequence == 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG':
            #    new_ubi_df = pd.read_csv(local_path + ".data/mol2_files/monomers/single_ubiquitin.txt", names=['#'])
            #elif next_working_FASTA_sequence == 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH':
            #    new_ubi_df = pd.read_csv(local_path + ".data/mol2_files/monomers/histag_ubiquitin.txt", names=['#'])



            ### download the correct ubiquiitin molecule
            #new_ubiAllAtoms, new_ubiAllBonds = get_BONDS_ATOMS(new_ubi_df)
            ### rotate new ubiquitin molecule
            #newUbiAtoms, conjugatingAtoms, newUbiBonds, conjugatingUbiBonds = rotation_transformation_of_new_ubiquitin(working_ubiAllAtoms, working_ubiAllBonds, new_ubiAllAtoms, new_ubiAllBonds, ubiquitinConjugationSite, conjugatingAtomsStart, conjugatingAtomsEnd)
            ### renumber the bonding for the ubiquitin molecule
            #concatUbiAtoms, concatUbiBonds = numbering_for_bonds_and_atoms(newUbiAtoms, conjugatingAtoms, newUbiBonds, conjugatingUbiBonds, ubiquitinConjugationSite, lastBackBoneSubstID)
            #


            ## ====================
            ## new ubiquitin addition in JSON, and in 
            #bra['children'], working_ubiAllAtoms, working_ubiAllBonds = inner_wrapper_relabelling_ubiquitin_numbers_and_JSON_to_MOL2(bra['children'], concatUbiAtoms, concatUbiBonds) 
            if first_dash:
                    first_dash = False
                    multimer_string_name = multimer_string_name + str(bra['site_name']) + '_'
            else: 
                multimer_string_name = multimer_string_name + '-' + str(bra['site_name']) + '_'
            bra['children'], multimer_string_name = inner_wrapper_getting_multimer_string_name(bra['children'], multimer_string_name) 
            multimer_string_name = multimer_string_name   

        logging.info(' ===== END OF LYSINE SITE =====  ')

        logging.info(' ===== END OF LYSINE SITE =====  ')

    multimer_string_name = multimer_string_name + ']'   
    logging.info(' ===== END OF PROTEIN - UBI NUMBER: ' + str(working_dictionary['chain_number']) + ' =====  ')

    return working_dictionary, multimer_string_name

##### functions to be applied to dataframe 
# global_deprot
# find_max_chain_number
def find_free_lysines(parent_dictionary):
    global chain_number_list
    global chain_length_list
    global free_lysine_list
    chain_number_list = [1]
    chain_length_list = []
    free_lysine_list = []

    if isinstance(parent_dictionary, str):
        x = parent_dictionary.replace("'", "\"")
        parent_dictionary = json.loads(x)

    logging.info(chain_number_list)
    adapted_dictionary = inner_wrapper_find_free_lysines(parent_dictionary)
    return adapted_dictionary, free_lysine_list
    #if not(specific_ubi_num in chain_number_list[:-1]): $
    #    raise TypeError('the ubiquitin specific is outside the range of the ubiquitin multimer')

### go into the protein, then the branching_sites... then loop the problem... 
### this function loops through the ubiquitin labeling the proteins...
def inner_wrapper_find_free_lysines(input_dictionary):
    ### super important because the dictionary is a mutable object
    working_dictionary = copy.deepcopy(input_dictionary)

    ## set the chain_number from the chain_number_list
    working_dictionary['chain_number'] = chain_number_list[-1]
    
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
            ### the difference.... 
            bra['children'] = inner_wrapper_find_free_lysines(bra['children']) 

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

def find_max_chain_number(parent_dictionary):
    global chain_number_list
    global chain_length_list
    global free_lysine_list
    chain_number_list = [1]
    chain_length_list = []
    free_lysine_list = []

    if isinstance(parent_dictionary, str):
        x = parent_dictionary.replace("'", "\"")
        parent_dictionary = json.loads(x)

    logging.info(chain_number_list)
    adapted_dictionary = inner_wrapper_find_max_chain_number(parent_dictionary)
    max_chain_number = chain_number_list[-1]-1
    return max_chain_number
    #if not(specific_ubi_num in chain_number_list[:-1]): $
    #    raise TypeError('the ubiquitin specific is outside the range of the ubiquitin multimer')

### go into the protein, then the branching_sites... then loop the problem... 
### this function loops through the ubiquitin labeling the proteins...
def inner_wrapper_find_max_chain_number(input_dictionary):
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
            ### the difference.... 
            bra['children'] = inner_wrapper_find_max_chain_number(bra['children']) 

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


def find_number_of_ABOC_SMAC(parent_dictionary):
    global chain_number_list
    global chain_length_list
    global free_lysine_list
    global ABOC_list
    global SMAC_list

    if isinstance(parent_dictionary, str):
        x = parent_dictionary.replace("'", "\"")
        parent_dictionary = json.loads(x)

    chain_number_list = [1]
    chain_length_list = []
    free_lysine_list = []
    ABOC_list = []
    SMAC_list = []
    logging.info(chain_number_list)
    adapted_dictionary = inner_wrapper_find_number_of_ABOC_SMAC(parent_dictionary)
    return ABOC_list, SMAC_list
    #if not(specific_ubi_num in chain_number_list[:-1]): $
    #    raise TypeError('the ubiquitin specific is outside the range of the ubiquitin multimer')

### go into the protein, then the branching_sites... then loop the problem... 
### this function loops through the ubiquitin labeling the proteins...
def inner_wrapper_find_number_of_ABOC_SMAC(input_dictionary):
    ### super important because the dictionary is a mutable object
    working_dictionary = copy.deepcopy(input_dictionary)

    ## set the chain_number from the chain_number_list
    working_dictionary['chain_number'] = chain_number_list[-1]
    
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
            if bra['children'] =='SMAC':
                SMAC_list.append([working_dictionary['chain_number'], str(bra['site_name']), str(bra['children'])])            
            if bra['children'] =='ABOC':
                ABOC_list.append([working_dictionary['chain_number'], str(bra['site_name']), str(bra['children'])])


        
        ### if the site is a protein.... 
        elif (bra['children'] == "") & (bra['site_name'] in ['K11','K6','K27','K29','K33','M1']): 
            logging.info("There is no Protecting Group on: " + str(bra['site_name']))

        elif (bra['children'] == "") & (bra['site_name'] in ['K48', 'K63']): 
            free_lysine_list.append([working_dictionary['chain_number'], str(bra['site_name'])])
            logging.info("There is no Protecting Group on: " + str(bra['site_name']))    
        
        elif isinstance(bra['children'], dict): 
            ### recursive function that calls itself... 
            logging.info('NEXT CHAIN: ' + str(bra['children']))
            ### the difference.... 
            bra['children'] = inner_wrapper_find_number_of_ABOC_SMAC(bra['children']) 

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

def find_number_of_ABOC(parent_dictionary):
    return len(find_number_of_ABOC_SMAC(parent_dictionary)[0])

def find_number_of_SMAC(parent_dictionary):
    return len(find_number_of_ABOC_SMAC(parent_dictionary)[1])



### you'll want to change this so that it is not recursive but pops up in each recursive loop
def validate_all_branching_sites(ubiquitin_structure):
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
            f"❌ Missing sites {missing_sites} in Ubiquitin {chain_number}."
        )

        # Recursively validate nested children
        for site in ubiquitin_dict["branching_sites"]:
            if isinstance(site["children"], dict):
                _check_sites(site["children"])

    _check_sites(ubiquitin_structure)