import logging

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

