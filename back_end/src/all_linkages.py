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
from src.main import * 
from tests.test_data import ubi_ubq_1

def process_branch_all(branch, working_dictionary, context):
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
    elif (branch['children'] == "") & (branch['site_name'] in ['M1']): 
        logging.info(f"Free lysine {branch['site_name']} found in chain {working_dictionary['chain_number']}.")
        
    # Handle 'K63', 'K48', 'K33', 'K29', 'K27','K11', 'K6'
    elif (branch['children'] == "") & (branch['site_name'] in ['K63', 'K48', 'K33', 'K29', 'K27','K11', 'K6']): 
        context["free_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name'])]]

    # Handle branches that have proteins bound  
    elif isinstance(branch["children"], dict):
        context["multimer_string_name"] += f"<{branch['site_name']}_"
        context["conjugated_lysines"] += [[working_dictionary['chain_number'], str(branch['site_name']), chain_number_list[-1]]]
        branch["children"], context = inner_wrapper_iterate_through_ubiquitin_all(
            branch["children"], context
        )
        context["multimer_string_name"] += ">"

    return branch, working_dictionary, context


def iterate_through_ubiquitin_all(parent_dictionary):
    
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

    adapted_dictionary, context = inner_wrapper_iterate_through_ubiquitin_all(
        parent_dictionary, context
    )

    return adapted_dictionary, context

def inner_wrapper_iterate_through_ubiquitin_all(input_dictionary, context):
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
        branch, working_dictionary, context = process_branch_all(branch, working_dictionary, context)

        # log end of branching details
        log_end_of_branching()

    log_end_of_protein(working_dictionary)

    # End of multimer string editing 
    context["multimer_string_name"] += ")"

    # Find the maximum chain number
    context = add_max_chain_number(context)
    
    return working_dictionary, context


def ubiquitin_building_all(
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
        "ubiquitin_added": False,  # Track if a ubiquitin has been added
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
        
    def inner_wrapper_ubiquitin_building_all(
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

            if (bra['site_name'] == lysine_residue) & (working_dictionary['chain_number'] == ubiquitin_number) & (context['ubiquitin_added'] == False):
                # Apply modification logic to the lysine site
                bra, working_dictionary = handle_lysine_modification(
                    bra, working_dictionary, ubi_molecule_to_add, ubiquitin_number, lysine_residue
                )
                context['ubiquitin_added'] = True  # Mark that a ubiquitin has been added so no more are further added

            # Log the state of the current branch
            if bra['children'] in ('SMAC', 'ABOC'):
                logging.info(f"Protecting Group: {bra['children']}")
            elif bra['children'] == "":
                logging.info(f"There is no Protecting Group on: {bra['site_name']}")
            elif isinstance(bra['children'], dict):
                logging.info(f"NEXT CHAIN: {bra['children']}")
                # Recursively process the next ubiquitin chain
                bra['children'], context = inner_wrapper_ubiquitin_building_all(
                    bra['children'], ubi_molecule_to_add, ubiquitin_number, lysine_residue, context
                )
            log_end_of_branching()

        log_end_of_protein(working_dictionary)

        return working_dictionary, context
    
    # Build structure recursively
    output_dictionary, context = inner_wrapper_ubiquitin_building_all(
        parent_dictionary, ubi_molecule_to_add, ubiquitin_number, lysine_residue, context
    )
    
    output_dictionary, output_context = iterate_through_ubiquitin_all(output_dictionary)

    return output_dictionary, output_context

def initialize_multimer_dicts_all(initial_acceptor):
    """
    Set up the initial multimer and context dictionaries.
    """
    multimer_dicts = {
        'multimers': [],
        'contexts': []
    }

    # Initialize acceptor and context
    acceptor, context = iterate_through_ubiquitin_all(initial_acceptor)

    # Add to the dictionary
    multimer_dicts['multimers'].append(acceptor)
    multimer_dicts['contexts'].append(context)

    return multimer_dicts


def defining_json_multimers_all(multimer_dicts, unprotected_ubi):
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
            working_multimer, working_context = ubiquitin_building_all(
                working_multimer,
                copy.deepcopy(unprotected_ubi),
                free_lysine[0],
                free_lysine[1]
            )

            new_multimer_dicts['multimers'].append(working_multimer)
            new_multimer_dicts['contexts'].append(working_context)

    return new_multimer_dicts

# ===========================================
# Functions for checking isomorphism in graphs
# For example, counting n-level topologies in higher-level graphs
# ===========================================
import networkx as nx
from networkx.algorithms import isomorphism
from itertools import combinations

def load_multimer_contexts(project_root: Path, multimer_size: int) -> dict:
    """
    Loads the contents of X_multimers_contexts.json from the specified project root,
    where X is the multimer size (e.g., 2 for dimers, 3 for trimers, etc.).
    
    Parameters:
        project_root (Path): Path object pointing to the project root directory.
        multimer_size (int): Size of the multimer (used to build the filename).
    
    Returns:
        dict: Parsed JSON data from the file.
    """
    filename = f"{multimer_size}_multimers_contexts.json"
    json_path = project_root / 'back_end' / 'data' / 'all_jsons' / filename

    with open(json_path, 'r') as f:
        data = json.load(f)
    
    return data

def n_in_higher_level(higher_level_edges, n_level_edges):
    """
    Counts the number of subgraph isomorphisms where a n-level topology
    (4-node graph) appears in a larger graphs.
    Avoids counting duplicate isomorphic matches by using frozenset node identifiers.
    """

    n_level_size = len(n_level_edges) + 1  # Assuming edges are in the form (source, label, target)

    def build_graph(edges):
        """
        Build a directed graph from a list of edges.
        Each edge is a tuple/list: (source, label, target).
        """
        G = nx.DiGraph()
        for src, label, dst in edges:
            G.add_edge(src, dst, label=label)
        return G

    G_pentamer = build_graph(higher_level_edges)
    G_tetramer = build_graph(n_level_edges)

    match_count = 0
    seen_matches = []

    for node_combo in combinations(G_pentamer.nodes, n_level_size):
        subgraph = G_pentamer.subgraph(node_combo).copy()

        matcher = isomorphism.DiGraphMatcher(
            subgraph,
            G_tetramer,
            edge_match=isomorphism.categorical_edge_match("label", None)
        )

        if matcher.is_isomorphic():
            node_set = frozenset(matcher.mapping.keys())
            if node_set not in seen_matches:
                seen_matches.append(node_set)
                match_count += 1

    return match_count

def get_multimer_edges_by_lysines(context_data: dict, 
                                  lysine_ids: dict
                                  ) -> dict:
    """
    Extracts edges from the context data where the source or target lysine ID is in the specified set.
    
    Parameters:
        context_data (dict): Dictionary containing multimer contexts.
        lysine_ids (set): Set of lysine IDs to filter edges by.
    
    Returns:
        dict: Dictionary with keys as multimer IDs and values as lists of edges matching the lysine IDs.
    """

    # DONE reveal all the edges for trimers 
    def reveal_edges(context_):
        """
        Gives the ID of the linkages from the context.
        """
        edges = context_['conjugated_lysines']
        return edges

    reveal_all_edges = {key: reveal_edges(value) for key, value in context_data.items()}

    # DONE only select edges with K63 and K48 linkages
    # Keep only entries where all edge labels are either K63 or K48
    def filter_edges_by_labels(reveal_all_edges, labels=lysine_ids):
        """
        Filters the edges in the reveal_all_edges dictionary to keep only those
        where all edge labels are in the specified set of labels.
        
        Parameters:
            reveal_all_edges (dict): Dictionary containing edges with their labels.
            labels (set): Set of labels to filter by (default is {'K63', 'K48'}).
        
        Returns:
            dict: Filtered dictionary with edges that match the criteria.
        """
        return {
            key: edges
            for key, edges in reveal_all_edges.items()
            if all(edge[1] in labels for edge in edges)
        } 

    return filter_edges_by_labels(reveal_all_edges)

def analyze_subgraph_containment(higher_level_dict, n_level_dict, progress_callback=None):
    """
    Analyzes the containment of n-level topologies within higher-level graphs.
    Counts how many times each n-level topology appears in each higher-level graph.
    Parameters:
        higher_level_dict (dict): Dictionary containing higher-level graph edges.
        n_level_dict (dict): Dictionary containing n-level graph edges.
        progress_callback (callable): Optional callback function to report progress.
    Returns:
        dict: A dictionary where keys are higher-level graph representations
              and values are dictionaries with n-level graph representations
              and their respective counts of occurrences.
    """
    results = {}
    import time
    
    total_items = len(higher_level_dict)
    sample_size = min(10, total_items)
    
    # Time the first 10 iterations for estimation
    start_time = time.time()

    for index, (high_key, high_edges) in enumerate(higher_level_dict.items()):
        high_str = str(high_edges)
        results[high_str] = {}
        for n_key, n_edges in n_level_dict.items():
            n_str = str(n_edges)
            count = n_in_higher_level(high_edges, n_edges)
            results[high_str][n_str] = count
        
        # Send progress update via callback
        if progress_callback:
            progress_callback({
                "type": "progress",
                "current": index + 1,
                "total": total_items,
                "message": f"Processed {high_str} with {len(n_level_dict)} n-level structures."
            })
        
        # After first 10 iterations, estimate total time
        if index + 1 == sample_size:
            elapsed_time = time.time() - start_time
            avg_time_per_iteration = elapsed_time / sample_size
            estimated_total_time = avg_time_per_iteration * total_items
            remaining_time = estimated_total_time - elapsed_time
            
            timing_info = {
                "type": "timing_analysis",
                "completed_iterations": sample_size,
                "elapsed_time": elapsed_time,
                "avg_time_per_iteration": avg_time_per_iteration,
                "estimated_total_time": estimated_total_time,
                "estimated_remaining_time": remaining_time,
                "estimated_total_seconds": estimated_total_time,
                "estimated_remaining_seconds": remaining_time
            }
            
            if progress_callback:
                progress_callback(timing_info)

    # Send completion notification
    if progress_callback:
        progress_callback({
            "type": "complete",
            "message": "Analysis completed successfully",
            "total_results": len(results)
        })
    
    return results

# =========================================================
# Build base multimers of polyubiquitin with all linkage
# This section initializes the multimers and builds them up to size 6.
# =========================================================

def build_all_linkages_multimers(
        monomer = ubi_ubq_1,
        largest_multimer_size=5,
        project_root_=project_root
        ):
    
    """Builds multimers of polyubiquitin with all linkages. 
    Args:
        ubi_ubq_1: The base ubiquitin structure.
        largest_multimer_size: The maximum size of the multimers to build.
        project_root_: The root path of the project for saving output files.
    """
    
    # Initialize the multimers
    multimers = initialize_multimer_dicts_all(monomer)

    # Build multimers of increasing size
    for multimer_size in range(2, largest_multimer_size+1):
        # Expand the multimer list by adding new ubiquitins
        multimers = defining_json_multimers_all(multimers, monomer)

        print(f"Multimer size {multimer_size} built with {len(multimers['multimers'])} entries.")

        # Remove duplicate entries from the multimer dictionary
        delete_duplicate_multimers(multimers)

        print(f"Following deletion, multimer size {multimer_size} built with {len(multimers['multimers'])} entries.")

        # Create output directory for JSON files
        output_dir = (
            project_root_
            / 'back_end'
            / 'data'
            / 'all_jsons'
        )
        output_dir.mkdir(parents=True, exist_ok=True)

        import json
        
        # Separate contexts and jsons from multimers
        contexts_data = {}
        jsons_data = {}
        
        if isinstance(multimers, dict):
            for key, value in multimers.items():
                if isinstance(value, list) and len(value) > 0:
                    # Check if this looks like context data (nested protein structures)
                    if isinstance(value[0], dict) and 'protein' in value[0] and 'branching_sites' in value[0]:
                        jsons_data[key] = value
                    # Check if this looks like json data (multimer information)
                    elif isinstance(value[0], dict) and 'chain_number_list' in value[0]:
                        contexts_data[key] = value
                    else:
                        # Default to jsons if unclear
                        jsons_data[key] = value
                else:
                    jsons_data[key] = value
        
        # Number the files as 2, 3, 4, 5 for multimer sizes 2, 3, 4, 5
        json_number = multimer_size
        
        # Save contexts data
        if contexts_data:
            contexts_numbered = {}
            counter = 1
            for key, value_list in contexts_data.items():
                if isinstance(value_list, list):
                    for individual_value in value_list:
                        contexts_numbered[str(counter)] = individual_value
                        counter += 1
                else:
                    contexts_numbered[str(counter)] = value_list
                    counter += 1
            
            contexts_path = output_dir / f"{json_number}_multimers_contexts.json"
            with open(contexts_path, 'w') as json_file:
                json.dump(contexts_numbered, json_file, indent=2)

        # Save jsons data
        if jsons_data:
            jsons_numbered = {}
            counter = 1
            for key, value_list in jsons_data.items():
                if isinstance(value_list, list):
                    for individual_value in value_list:
                        jsons_numbered[str(counter)] = individual_value
                        counter += 1
                else:
                    jsons_numbered[str(counter)] = value_list
                    counter += 1
            
            jsons_path = output_dir / f"{json_number}_multimers_jsons.json"
            with open(jsons_path, 'w') as json_file:
                json.dump(jsons_numbered, json_file, indent=2)

# =========================================================