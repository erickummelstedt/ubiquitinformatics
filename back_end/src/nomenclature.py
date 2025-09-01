from pathlib import Path
import pandas as pd
import sys
import copy

# Dynamically get the backend path relative to this file
current_file = Path(__file__).resolve()
project_root = current_file.parents[2]  # Go up to project root
sys.path.insert(0, str(project_root))
local_path = project_root / 'back_end'
sys.path.insert(0, str(local_path))

from src.main import *
from src.simulation import *
from src.utils.utils import *
from src.data_cleaning import *
from src.all_linkages import *
from tests.test_data import *

"""
UBIQUITIN NOMENCLATURE SYSTEM
============================

This module implements a comprehensive nomenclature system for representing ubiquitin chain structures 
in multiple complementary formats. The system provides bidirectional conversion between different 
notation styles to support various scientific and computational use cases.

NOMENCLATURE FORMATS:
--------------------

1. TREE NOMENCLATURE (A1B2C4D6 format):
   - Hierarchical representation using letters (A,B,C,D...) for levels and numbers for positions
   - Level System: A = level 1, B = level 2, C = level 3, D = level 4, etc.
   - Position System: Each level has 2^(level-1) capacity (A:1, B:2, C:4, D:8...)
   - Linkage Encoding by Position Number:
     * Position 1 = K63 linkage
     * Position 2 = K48 linkage  
     * Position 3 = K33 linkage
     * Position 4 = K29 linkage
     * Position 5 = K27 linkage
     * Position 6 = K11 linkage
     * Position 7 = K6 linkage
     * Position 8 = K63 linkage (cycle repeats)
     * Pattern: [K63, K48, K33, K29, K27, K11, K6] repeats for positions 1-7, 8-14, 15-21, etc.
   - Binding Rules: Position x can bind to positions (7x-6) through (7x) in next level
     * Example: A1 → B1,B2,B3,B4,B5,B6,B7 (positions 1-7)
     * Example: A2 → B8,B9,B10,B11,B12,B13,B14 (positions 8-14)
   - Example: A1B1B2C1C4 represents a branched tetramer with specific connectivity

2. COMPACT EDGE FORMAT (1A2-2A3-3A4 format):
   - Linear representation of direct connections between ubiquitin units
   - Format: SourceUnitLetterDestinationUnit (e.g., 1A2 = unit 1 connects to unit 2 via K63)
   - Letter Mapping: A=K63, B=K48, C=K33, D=K29, E=K27, F=K11, G=K6
   - Separators: Supports both dash (-) and comma (,) separation
   - Example: "1A2-2A3-3B4" = linear chain with K63 linkages and one K48 branch

3. CONJUGATED LYSINES FORMAT ([[1,'K63',2],...] format):
   - Explicit list representation showing source, linkage type, and destination
   - Format: [[source_unit, 'Kxx', destination_unit], ...]
   - Preserves preorder traversal information for tree reconstruction
   - Example: [[1,'K63',2], [2,'K63',3], [1,'K48',4]] = branched structure

4. FORMATTED EDGES (1 → K63 → 2 format):
   - Human-readable arrow notation for visualization
   - Shows explicit linkage sites and connection directions
   - Used primarily for display and user interface components
   - Example: "1 → K63 → 2, 2 → K63 → 3, 1 → K48 → 4"

CONVERSION CAPABILITIES:
-----------------------
- Tree Nomenclature ↔ Conjugated Lysines
- Tree Nomenclature → Polyubiquitin Structure (via ubiquitin_building_all)
- Compact Edges → Conjugated Lysines → Tree Nomenclature
- All formats → Formatted Edges for display

VALIDATION FEATURES:
-------------------
- Tree structure validation (binding rules, level constraints)
- Preorder traversal preservation during conversions
- Error handling for malformed inputs
- Consistent position calculation algorithms

SCIENTIFIC APPLICATIONS:
-----------------------
- Automated synthesis planning for complex ubiquitin chains
- Structure-function relationship analysis
- Computational modeling of ubiquitin signaling pathways
- Standardized notation for research documentation and databases
"""

### BETTER VERSION ### 

def format_edges(edges):
    """
    Convert edges with lysine labels into compact string form using fixed mapping:
        K63 -> A, K48 -> B, K33 -> C, K29 -> D, K27 -> E, K11 -> F, K6 -> G
    Example:
        [[1, 'K63', 2], [1, 'K48', 4]] -> "1A2, 1B4"
    """
    fixed_map = {
        'K63': 'A',
        'K48': 'B',
        'K33': 'C',
        'K29': 'D',
        'K27': 'E',
        'K11': 'F',
        'K6':  'G'
    }
    return '-'.join(f"{s}{fixed_map.get(k, '?')}{d}" for s, k, d in edges)


import re

# Fixed, canonical mapping
LETTER_TO_LYS = {
    'A': 'K63',
    'B': 'K48',
    'C': 'K33',
    'D': 'K29',
    'E': 'K27',
    'F': 'K11',
    'G': 'K6',
}


def parse_compact_edges(compact):
    """
    Parse compact edge notation into [[src, 'Kxx', dst], ...].

    Accepts either:
      - a single string like "1A2, 1B4, 2A3, 4B5" (comma-separated)
      - a single string like "1A2-2A3-3A4-4B5" (dash-separated)
      - an iterable of strings like ["1A2", "1B4", "2A3", "4B5"]

    Returns: list of [int, str, int], e.g. [[1,'K63',2], ...]
    """
    # Normalize to a list of tokens like ["1A2", "1B4", ...]
    if isinstance(compact, str):
        # split on commas or dashes; also allow arbitrary whitespace
        if ',' in compact:
            tokens = [t.strip() for t in compact.split(',') if t.strip()]
        elif '-' in compact:
            tokens = [t.strip() for t in compact.split('-') if t.strip()]
        else:
            # Single token
            tokens = [compact.strip()] if compact.strip() else []
    else:
        tokens = [str(t).strip() for t in compact if str(t).strip()]

    edges = []
    pat = re.compile(r'^(\d+)\s*([A-G])\s*(\d+)$')
    for tok in tokens:
        m = pat.match(tok)
        if not m:
            raise ValueError(f"Invalid compact edge token: {tok!r} (expected like '12A7')")
        s, letter, d = m.groups()
        lys = LETTER_TO_LYS.get(letter)
        if lys is None:
            raise ValueError(f"Unknown lysine letter: {letter!r}")
        edges.append([int(s), lys, int(d)])
    return edges

# NEW FUNCTION: build_polyubiquitin_from_edges
def build_polyubiquitin_from_edges(connections):
    """
    Build a polyubiquitin structure from an explicit edge list.

    Args:
        connections (list): List like [[src, 'Kxx', dst], ...],
                            e.g. [[1, 'K63', 2], [1, 'K48', 4], [2, 'K63', 3], [4, 'K48', 5]]

    Returns:
        dict: The final polyubiquitin structure built from the edges.
    """
    print("hello")
    if not connections:
        # Nothing to add: just iterate the base structure
        output_structure, _ = iterate_through_ubiquitin(histag_ubi_ubq_1)
        return output_structure

    print("hello2")
    # Start with the base ubiquitin molecule
    current_structure = histag_ubi_ubq_1

    # Apply each connection iteratively; dst is implied/unused by ubiquitin_building_all
    for edge in connections:
        try:
            src, lys, dst = edge  
        except ValueError:
            # Skip malformed edges
            continue

        # Normalize types
        try:
            ubiquitin_number = int(src)
        except (TypeError, ValueError):
            # If src isn't an int, skip this edge
            continue
        lysine_residue = str(lys)

        print("hello3")

        # Prepare the new ubiquitin to add. give it the correct chain number
        new_ubi_ubq_1 = ubi_ubq_1
        new_ubi_ubq_1= convert_json_to_dict(new_ubi_ubq_1)
        new_ubi_ubq_1['chain_number'] = int(dst)

        # Build the chain
        current_structure, _ = ubiquitin_building_all(
            parent_dictionary=current_structure,
            ubi_molecule_to_add=new_ubi_ubq_1,
            ubiquitin_number=ubiquitin_number,
            lysine_residue=lysine_residue,
        )

    output_structure, _ = iterate_through_ubiquitin(current_structure)
    return output_structure


import re

def multimer_length_from_nomenclature(compact: str) -> int:
    """
    Return the number of unique ubiquitin units in a compact edge string.

    Example:
        "1A2-2A3-3A4-4B5" -> 5  (units: {1,2,3,4,5})

    Accepts letters A–G for lysines (fixed mapping), ignores whitespace.
    """
    # Find all "number Letter number" patterns, e.g. 1A2, 12B7, etc.
    pairs = re.findall(r'(\d+)\s*[A-G]\s*(\d+)', compact)
    nodes = set()
    for s, d in pairs:
        nodes.add(int(s))
        nodes.add(int(d))
    return len(nodes)



# ==================================
# ==== Jeff K48-K63 Nomenclature Conversion
# ==================================

def conjugated_lysines_to_jeff_K48_K63_nomenclature(conjugated_lysines):
    """
    Convert conjugated lysines format to tree nomenclature (K48-K63 ONLY)
    
    Args:
        conjugated_lysines (list): List of connections in format [[src, linkage, dst], ...]
                                  Example: [[1, 'K63', 2], [2, 'K63', 3], [1, 'K48', 4], [4, 'K48', 5]]
                                  Order matters - represents preorder traversal
                                  ONLY K48 and K63 linkages are supported
    
    Returns:
        str: Tree nomenclature string (e.g., 'A1B1B2C3D5')
    
    Raises:
        ValueError: If any linkage other than K48 or K63 is found
    """
    if not conjugated_lysines:
        return "A1"  # Single ubiquitin
    
    # Validate that only K48 and K63 linkages are present
    allowed_lysines = {'K48', 'K63'}
    for src, linkage, dst in conjugated_lysines:
        if linkage not in allowed_lysines:
            raise ValueError(f"K48-K63 nomenclature only supports K48 and K63 linkages. Found: {linkage}")
    
    # Build mapping from node number to the order it appears in connections
    node_to_order = {}
    all_nodes = set()
    
    # Track the order nodes appear in the conjugated lysines list
    order_counter = 0
    for src, linkage, dst in conjugated_lysines:
        if src not in node_to_order:
            node_to_order[src] = order_counter
            order_counter += 1
        if dst not in node_to_order:
            node_to_order[dst] = order_counter
            order_counter += 1
        all_nodes.add(src)
        all_nodes.add(dst)
    
    # Build tree structure
    children = {}  # parent -> list of (child, linkage, order)
    parents = {}   # child -> (parent, linkage)
    
    for src, linkage, dst in conjugated_lysines:
        if src not in children:
            children[src] = []
        children[src].append((dst, linkage, node_to_order[dst]))
        parents[dst] = (src, linkage)
    
    # Find root node (node with no parent)
    root = None
    for node in all_nodes:
        if node not in parents:
            root = node
            break
    
    if root is None:
        return "Error: No root node found"
    
    # Assign positions using the actual order from conjugated lysines
    nomenclature_map = {}  # node -> (level_letter, position)
    
    def assign_positions(node, level, position):
        level_letter = chr(ord('A') + level)
        nomenclature_map[node] = (level_letter, position)
        
        if node in children:
            # Sort children by their appearance order in the original list
            # This preserves the preorder traversal where K63 branches are visited before K48
            sorted_children = sorted(children[node], key=lambda x: x[2])  # Sort by order
            
            for child, linkage, _ in sorted_children:
                # Calculate child position based on parent position and linkage type
                # ONLY K63 and K48 are allowed here due to validation above
                if linkage == 'K63':
                    child_position = 2 * position - 1
                elif linkage == 'K48':
                    child_position = 2 * position
                else:
                    # This should never happen due to validation, but just in case
                    raise ValueError(f"Unexpected linkage type in K48-K63 nomenclature: {linkage}")
                
                assign_positions(child, level + 1, child_position)
    
    # Start DFS from root with position 1
    assign_positions(root, 0, 1)
    
    # Build nomenclature string by collecting all positions and sorting by level then position
    all_positions = []
    for node in nomenclature_map:
        level_letter, position = nomenclature_map[node]
        level_num = ord(level_letter) - ord('A')
        all_positions.append((level_num, position, level_letter, node))
    
    # Sort by level first, then by position
    all_positions.sort(key=lambda x: (x[0], x[1]))
    
    # Build the nomenclature string
    nomenclature_parts = []
    for level_num, position, level_letter, node in all_positions:
        nomenclature_parts.append(f"{level_letter}{position}")
    
    return ''.join(nomenclature_parts)


# ==================================
# ==== Jeff all lysines Nomenclature Conversion
# ==================================

def conjugated_lysines_to_jeff_all_lysines_nomenclature(conjugated_lysines):
    """
    Convert conjugated lysines format to tree nomenclature using 7-lysine mapping system
    
    Args:
        conjugated_lysines (list): List of connections in format [[src, linkage, dst], ...]
                                  Example: [[1, 'K63', 2], [2, 'K33', 3], [1, 'K48', 4]]
                                  Order matters - represents preorder traversal
    
    Returns:
        str: Tree nomenclature string using 7-lysine cyclic position mapping
    
    Position Mapping:
        K63: positions 1, 8, 15, 22... (1 + 7*n)
        K48: positions 2, 9, 16, 23... (2 + 7*n)  
        K33: positions 3, 10, 17, 24... (3 + 7*n)
        K29: positions 4, 11, 18, 25... (4 + 7*n)
        K27: positions 5, 12, 19, 26... (5 + 7*n)
        K11: positions 6, 13, 20, 27... (6 + 7*n)
        K6:  positions 7, 14, 21, 28... (7 + 7*n)
    
    Binding Rules: Parent at position x can bind to positions (7*x-6) through (7*x) in next level
    """
    if not conjugated_lysines:
        return "A1"  # Single ubiquitin
    
    # Lysine to base position mapping (1-7 cycle)
    lysine_to_base_position = {
        'K63': 1,
        'K48': 2,
        'K33': 3,
        'K29': 4,
        'K27': 5,
        'K11': 6,
        'K6': 7
    }
    
    # Build mapping from node number to the order it appears in connections
    node_to_order = {}
    all_nodes = set()
    
    # Track the order nodes appear in the conjugated lysines list
    order_counter = 0
    for src, linkage, dst in conjugated_lysines:
        if src not in node_to_order:
            node_to_order[src] = order_counter
            order_counter += 1
        if dst not in node_to_order:
            node_to_order[dst] = order_counter
            order_counter += 1
        all_nodes.add(src)
        all_nodes.add(dst)
    
    # Build tree structure
    children = {}  # parent -> list of (child, linkage, order)
    parents = {}   # child -> (parent, linkage)
    
    for src, linkage, dst in conjugated_lysines:
        if src not in children:
            children[src] = []
        children[src].append((dst, linkage, node_to_order[dst]))
        parents[dst] = (src, linkage)
    
    # Find root node (node with no parent)
    root = None
    for node in all_nodes:
        if node not in parents:
            root = node
            break
    
    if root is None:
        return "Error: No root node found"
    
    # Assign positions using the 7-lysine mapping system
    nomenclature_map = {}  # node -> (level_letter, position)
    
    def assign_positions(node, level, position):
        level_letter = chr(ord('A') + level)
        nomenclature_map[node] = (level_letter, position)
        
        if node in children:
            # Sort children by their appearance order in the original list
            # This preserves the preorder traversal
            sorted_children = sorted(children[node], key=lambda x: x[2])  # Sort by order
            
            # Calculate the range of positions this parent can bind to: (7*position-6) to (7*position)
            min_child_position = 7 * position - 6
            max_child_position = 7 * position
            
            # Track which positions are used to avoid conflicts
            used_positions = set()
            
            for child, linkage, _ in sorted_children:
                # Get base position for this lysine type (1-7)
                base_pos = lysine_to_base_position.get(linkage)
                if base_pos is None:
                    continue  # Skip unknown lysine types
                
                # Find the correct position within the allowed range
                # Start with base position and increment by 7 until we're in range
                child_position = base_pos
                
                # Adjust to be within the parent's binding range
                while child_position < min_child_position:
                    child_position += 7
                
                # If this position is beyond the range, use the highest available position
                if child_position > max_child_position:
                    # Find the largest valid position for this lysine type within range
                    for test_pos in range(max_child_position, min_child_position - 1, -1):
                        if (test_pos - base_pos) % 7 == 0 and test_pos not in used_positions:
                            child_position = test_pos
                            break
                    else:
                        # If no valid position found, skip this child
                        continue
                
                # Ensure position is not already used
                while child_position in used_positions and child_position <= max_child_position:
                    child_position += 7
                
                if child_position <= max_child_position:
                    used_positions.add(child_position)
                    assign_positions(child, level + 1, child_position)
    
    # Start DFS from root with position 1
    assign_positions(root, 0, 1)
    
    # Build nomenclature string by collecting all positions and sorting by level then position
    all_positions = []
    for node in nomenclature_map:
        level_letter, position = nomenclature_map[node]
        level_num = ord(level_letter) - ord('A')
        all_positions.append((level_num, position, level_letter, node))
    
    # Sort by level first, then by position
    all_positions.sort(key=lambda x: (x[0], x[1]))
    
    # Build the nomenclature string
    nomenclature_parts = []
    for level_num, position, level_letter, node in all_positions:
        nomenclature_parts.append(f"{level_letter}{position}")
    
    return ''.join(nomenclature_parts)

# ==================================
# ==== Numerical System Conversion
# ==================================

def tree_nomenclature_to_numerical_system(tree_nomenclature):
    """
    Convert tree nomenclature to a purely numerical system using base-7 hierarchy.
    
    Args:
        tree_nomenclature (str): Tree nomenclature string (e.g., 'A1B1B2C1C4')
    
    Returns:
        str: Comma-separated string of numerical values where each position is converted using:
             A = 0, B = 7, C = 7*7 = 49, D = 7*7*7 = 343, etc.
             Final value = level_multiplier + position_number
             Example: D9 = 7*7*7 + 9 = 343 + 9 = 352
    
    Example:
        'A1B1B2C1C4' -> "1, 8, 9, 50, 53"
        (A1=0+1=1, B1=7+1=8, B2=7+2=9, C1=49+1=50, C4=49+4=53)
    """
    if not tree_nomenclature:
        return ""
    
    import re
    
    # Find all level-position pairs (e.g., 'A1', 'B2', 'C4')
    # Only accept letters A-Z followed by digits
    pattern = r'([A-Z])(\d+)'
    matches = re.findall(pattern, tree_nomenclature)
    
    if not matches:
        return ""
    
    # Validate that the entire string consists only of valid nomenclature parts
    # Reconstruct what should be the input from the matches
    reconstructed = ''.join(f"{letter}{number}" for letter, number in matches)
    if reconstructed != tree_nomenclature:
        return ""  # Invalid format - contains extra characters
    
    numerical_values = []
    
    for level_letter, position_str in matches:
        # Calculate level multiplier: A=0, B=7^1, C=7^2, D=7^3, etc.
        level_index = ord(level_letter) - ord('A')
        level_multiplier = 7 ** level_index if level_index > 0 else 0
        
        # Get position number
        position_number = int(position_str)
        
        # Calculate final numerical value
        numerical_value = level_multiplier + position_number
        numerical_values.append(numerical_value)
    
    return ", ".join(map(str, numerical_values))



# =========================================================
# Mass Spec Dictionary Functions
# Extract unique FASTA sequences for mass spectrometry analysis
# =========================================================

def extract_fasta_sequences_for_mass_spec(parent_dictionary):
    """
    Traverse the polyubiquitin structure and extract unique FASTA sequences
    to build a mass spec dictionary mapping sequence names to their FASTA sequences.
    
    Args:
        parent_dictionary (dict or str): Ubiquitin structure as a dictionary or JSON string.
    
    Returns:
        dict: Dictionary mapping sequence names (e.g., "his-Ub", "Ub") to their FASTA sequences.
    """
    
    # Ensure that the parent dictionary is a JSON
    parent_dictionary = convert_json_to_dict(parent_dictionary)
    
    # Ensure that the parent dictionary has all the valid keys 
    validate_protein_keys(parent_dictionary)
    
    # Initialize context object to track unique FASTA sequences
    fasta_context = {
        "unique_sequences": {},  # Dictionary to store unique sequence_name -> FASTA_sequence mappings
        "chain_number_list": [1]
    }
    
    # Start the recursive traversal
    _, fasta_context = inner_wrapper_extract_fasta_sequences(parent_dictionary, fasta_context)
    
    return fasta_context["unique_sequences"]


def inner_wrapper_extract_fasta_sequences(input_dictionary, fasta_context):
    """
    Recursively process nested dictionaries to extract unique FASTA sequences
    and map them to appropriate sequence names.
    """
    working_dictionary = copy.deepcopy(input_dictionary)
    
    # Ensure that the working dictionary is a JSON
    working_dictionary = convert_json_to_dict(working_dictionary)
    
    # Ensure that the working dictionary has all the valid keys 
    validate_protein_keys(working_dictionary)
    
    # Process current protein and assign chain number
    working_dictionary, fasta_context = process_current_protein_fasta(working_dictionary, fasta_context)
    
    # Extract FASTA sequence and determine sequence name
    fasta_sequence = working_dictionary["FASTA_sequence"]
    
    # Determine sequence name based on FASTA sequence
    if fasta_sequence == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH":
        if working_dictionary["chain_number"] == 1:
            sequence_name = "his-Ub"
        else:
            sequence_name = "his-Ub"  # Keep consistent naming
    elif fasta_sequence == "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG":
        if working_dictionary["chain_number"] == 1:
            sequence_name = "Ub"
        else:
            sequence_name = "Ub"  # Keep consistent naming
    else:
        # For any other FASTA sequences, use a generic naming pattern
        sequence_name = f"Unknown_Protein_{working_dictionary['protein']}"
    
    # Add to unique sequences dictionary
    fasta_context["unique_sequences"][sequence_name] = fasta_sequence
    
    # Log protein details for debugging
    logging.info(f"Processing protein chain {working_dictionary['chain_number']}: {sequence_name}")
    logging.info(f"FASTA sequence: {fasta_sequence[:50]}..." if len(fasta_sequence) > 50 else f"FASTA sequence: {fasta_sequence}")
    
    # Process branching sites recursively
    working_branching_sites = working_dictionary.get("branching_sites", [])
    
    for branch in working_branching_sites:
        # Log branching details
        logging.info(f"Processing branch: {branch['site_name']} with children: {type(branch['children'])}")
        
        # If branch has a protein child, recursively process it
        if isinstance(branch["children"], dict):
            branch["children"], fasta_context = inner_wrapper_extract_fasta_sequences(
                branch["children"], fasta_context
            )
    
    return working_dictionary, fasta_context


def process_current_protein_fasta(working_dictionary, fasta_context):
    """
    Process the current protein for FASTA extraction, similar to process_current_protein
    but focused on sequence extraction needs.
    """
    # Set the current chain number from context
    working_dictionary['chain_number'] = fasta_context['chain_number_list'][-1]
    
    # Set chain length
    working_dictionary['chain_length'] = len(working_dictionary['FASTA_sequence'])
    
    # Increment chain_number for future recursive calls
    fasta_context['chain_number_list'].append(fasta_context['chain_number_list'][-1] + 1)
    
    return working_dictionary, fasta_context


def build_mass_spec_dictionary(parent_dictionary):
    """
    Convenience function to build a mass spec dictionary from a polyubiquitin structure.
    
    Args:
        parent_dictionary (dict or str): Ubiquitin structure as a dictionary or JSON string.
    
    Returns:
        dict: Dictionary with sequence names as keys and FASTA sequences as values.
              Example: {"his-Ub": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
                       "Ub": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"}
    """
    return extract_fasta_sequences_for_mass_spec(parent_dictionary)

# =========================================================