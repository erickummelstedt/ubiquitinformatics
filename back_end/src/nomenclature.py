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
from src.building_blocks import *

"""
UBIQUITIN NOMENCLATURE SYSTEM
============================

This module implements a comprehensive nomenclature system for representing ubiquitin chain structures 
in multiple complementary formats. The system provides bidirectional conversion between different 
notation styles to support various scientific and computational use cases.


NOMENCLATURE FORMATS:

Note: The vertices (V) are numbered by pre-order traversal, enabling a standardized and biologically consistent indexing scheme. Numbering proceeds from the C-terminus to the N-terminus, following the lysine order: K63, K48, K33, K29, K27, K11, K6.
In pre-order traversal, each branch is fully explored before the numbering returns to the current chain, ensuring that all descendants of a branch are numbered before continuing along the parent chain.

1. UbID: PREORDER NODE LABELING (where nodes 1 = A, 2 = B format):
    - Assigns each node a unique label based on preorder tree traversal. (A, B, C, ...)
    - Labels each edge according to the lysine linkage type.
    - Format: Node labels (A, B, C, ...) represent node numbers in preorder traversal (A=1, B=2, C=3, etc.) 
    - Format: Edge labels the number between the letters indicates the edge or lysine linkage (63=K63, 48=K48, etc.).
    - Highlights the connectivity and linkage chemistry between nodes.
    - Example: "A48B-B63C" = node A (1) connects to node B (2) via K48, node B (2) connects to node C (3) via K63.
    - Further Example: A63B-B63C-C48D-B48E = Pentamer trimer with 2 K48 branches at position B(2) and C(3)

2. PREORDER NODE LABELING (where edges K63 = A, K48 = B... format):
    - Assigns each node a unique label based on preorder tree traversal.(1, 2, 3, ...)
    - Format: Node labels (1, 2, 3, ...) represent node numbers in preorder traversal (A=1, B=2, C=3, etc.)
    - Format: Edge labels (A, B, C, ...) correspond to the edge or lysine linkage (A=K63, B=K48, etc.).
    - Highlights the connectivity and linkage chemistry between nodes.
    - Example: "1B2-2A3" = node 1 connects to node 2 via B (K48), node 2 connects to node 3 via A (K63).
    - Further Example: 1A2-2A3-3B4-2B5 = Pentamer trimer with 2 K48 branches at position 2 and 3
    
3. GRAPH-BASED NOMENCLATURE WITH PREORDER NUMBERING (Bode/Majima/Kummelstedt):
    - Represents the chain as a nested structure, preserving the order in which nodes are added.
    - Format: Ub1,63(Ub2,63(Ub3,48(Ub4)),48(Ub5))
    - Useful for isomorphism checks and structural classification, focusing on connectivity while nodes are still numbered by preorder traversal.
    - Example: Ub1,63(Ub2,63(Ub3,48(Ub4)),48(Ub5))

4. CHEMISTRY-STYLE ALL NODE LABELED NOMENCLATURE (Bode):
    - Uses a new chemical notation conventions to label all nodes in the chain.
    - Format: Each node is assigned a label based on its chemical context and position.
    Level System:
        Level 1 = A, Level 2 = B, Level 3 = C, Level 4 = D, etc.

    Position Mapping:
        K63: evens with uppercase letter (e.g., B2, C2, D4)
        K48: odds with uppercase letter (e.g., B1, C1, D3)
        K33: evens with uppercase letter + ' (e.g., B'2, C'2, D'4)
        K29: odds with uppercase letter + ' (e.g., B'1, C'1, D'3)
        K11: evens with lowercase letter (e.g., b2, c2, d4)
        K6:  odds with lowercase letter (e.g., b1, c1, d3)
        K27: evens with lowercase letter + ' (e.g., b'2, c'2, d'4)
        M1:  odds with lowercase letter + ' (e.g., b'1, c'1, d'3)

    Formula: position = ((parent_letter_size + parent_number) - 1) * 2 + child_number

    Where:
        parent_letter_size: Based on parent's notation type
            - 0 for uppercase without prime (A1, B2)
            - 2 for uppercase with prime (A'1, B'2)
            - 4 for lowercase without prime (a1, b2)
            - 6 for lowercase with prime (a'1, b'2)
        parent_number: Numeric part of parent's notation
        child_number: +1 for K48/K29/K6/M1, +2 for K63/K33/K11/K27

5. GRAPH-BASED NOMENCLATURE WITHOUT PREORDER (Strieter/Shestoperova/Ivanov):
    - Represents the chain as a graph without enforcing preorder traversal or numbering.
    - Format: Nested structure without explicit node numbers, e.g., Ub,63(Ub,63(Ub,48(Ub)),48(Ub))
    - Useful for isomorphism checks and structural classification, focusing on connectivity rather than order.
    - Example: Ub,63(Ub,63(Ub,48(Ub)),48(Ub))


Other
6. CONJUGATED LYSINES FORMAT ([[1,'K63',2],...] format):
    - Explicit list representation showing source, linkage type, and destination
    - Format: [[source_unit, 'Kxx', destination_unit], ...]
    - Preserves preorder traversal information for tree reconstruction
    - Example: [[1,'K63',2], [2,'K63',3], [1,'K48',4]] = branched structure

7. FORMATTED EDGES (1 → K63 → 2 format):
    - Human-readable arrow notation for visualization
    - Shows explicit linkage sites and connection directions
    - Used primarily for display and user interface components
    - Example: "1 → K63 → 2, 2 → K63 → 3, 1 → K48 → 4"


SCIENTIFIC APPLICATIONS:
-----------------------
- Automated synthesis planning for complex ubiquitin chains
- Structure-function relationship analysis
- Computational modeling of ubiquitin signaling pathways
- Standardized notation for research documentation and databases
"""

# ==================================
# === Translating Preorder ABC nomenclature to explicit edges
# ==================================

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
      - a single string like "A63B, A48D, B63C, D48E" (comma-separated)
      - a single string like "A63B-A48D-B63C-D48E" (dash-separated)
      - an iterable of strings like ["A63B", "A48D", "B63C", "D48E"]

    Format: Letter-LysineNumber-Letter (e.g., A63B = node A connects to node B via K63)
    Node mapping: A=1, B=2, C=3, D=4, etc.
    Lysine mapping: 63=K63, 48=K48, 33=K33, 29=K29, 27=K27, 11=K11, 6=K6

    Returns: list of [int, str, int], e.g. [[1,'K63',2], ...]
    """
    # Normalize to a list of tokens like ["A63B", "A48D", ...]
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

    # Reverse mapping from lysine numbers to K notation
    LYSINE_NUM_TO_K = {
        '63': 'K63',
        '48': 'K48',
        '33': 'K33',
        '29': 'K29',
        '27': 'K27',
        '11': 'K11',
        '6': 'K6'
    }
    
    def letter_to_number(letter):
        """Convert letter to number: A->1, B->2, C->3, etc."""
        return ord(letter.upper()) - ord('A') + 1

    edges = []
    pat = re.compile(r'^([A-Z])(\d+)([A-Z])$')
    for tok in tokens:
        m = pat.match(tok)
        if not m:
            raise ValueError(f"Invalid compact edge token: {tok!r} (expected like 'A63B')")
        src_letter, lysine_num, dst_letter = m.groups()
        
        # Convert letters to numbers
        src_num = letter_to_number(src_letter)
        dst_num = letter_to_number(dst_letter)
        
        # Convert lysine number to K notation
        lys = LYSINE_NUM_TO_K.get(lysine_num)
        if lys is None:
            raise ValueError(f"Unknown lysine number: {lysine_num!r}")
        
        edges.append([src_num, lys, dst_num])
    return edges

# build_polyubiquitin_from_edges
def build_polyubiquitin_from_edges(connections):
    """
    Build a polyubiquitin structure from an explicit edge list.

    Args:
        connections (list): List like [[src, 'Kxx', dst], ...],
                            e.g. [[1, 'K63', 2], [1, 'K48', 4], [2, 'K63', 3], [4, 'K48', 5]]

    Returns:
        dict: The final polyubiquitin structure built from the edges.
    """

    # Start with the base ubiquitin molecule
    new_ubi_ubq_1 = ubi_ubq_1
    new_ubi_ubq_1= convert_json_to_dict(new_ubi_ubq_1)
    new_ubi_ubq_1['chain_number'] = int(1)
    current_structure = new_ubi_ubq_1.copy()

    if not connections:
        # Nothing to add: just iterate the base structure
        output_structure, _ = iterate_through_ubiquitin(ubi_ubq_1)
        return output_structure

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

# build_polyubiquitin_from_edges
def build_polyubiquitin_from_edges_with_histag(connections):
    """
    Build a polyubiquitin structure from an explicit edge list.

    Args:
        connections (list): List like [[src, 'Kxx', dst], ...],
                            e.g. [[1, 'K63', 2], [1, 'K48', 4], [2, 'K63', 3], [4, 'K48', 5]]

    Returns:
        dict: The final polyubiquitin structure built from the edges.
    """

    # Start with the base ubiquitin molecule
    new_ubi_ubq_1 = histag_ubi_ubq_1
    new_ubi_ubq_1= convert_json_to_dict(new_ubi_ubq_1)
    new_ubi_ubq_1['chain_number'] = int(1)
    current_structure = new_ubi_ubq_1.copy()

    if not connections:
        # Nothing to add: just iterate the base structure
        output_structure, _ = iterate_through_ubiquitin(ubi_ubq_1)
        return output_structure

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
        "A63B-A48C-B63D" -> 4  (units: {A,B,C,D} = {1,2,3,4})

    Accepts the A63B format (Letter-LysineNumber-Letter), ignores whitespace.
    """
    # Find all "Letter LysineNumber Letter" patterns, e.g. A63B, B48C, etc.
    pairs = re.findall(r'([A-Z])\d+([A-Z])', compact)
    nodes = set()
    
    def letter_to_number(letter):
        """Convert letter to number: A->1, B->2, C->3, etc."""
        return ord(letter.upper()) - ord('A') + 1
    
    for src_letter, dst_letter in pairs:
        nodes.add(letter_to_number(src_letter))
        nodes.add(letter_to_number(dst_letter))
    return len(nodes)


# ==================================
# ==== Nomenclature functions
# =================================

def format_nomenclature_preorder_1A2(edges):
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


def format_nomenclature_preorder_A63B(edges):
    """
    Convert edges with lysine labels into Jeff's preorder format using lysine numbers and node letters.
    Lysine mapping: K63 -> 63, K48 -> 48, K33 -> 33, K29 -> 29, K27 -> 27, K11 -> 11, K6 -> 6
    Node mapping: 1 -> A, 2 -> B, 3 -> C, 4 -> D, etc.
    
    Example:
        [[1, 'K63', 2], [1, 'K48', 3]] -> "A63B-A48C"
    """
    lysine_map = {
        'K63': '63',
        'K48': '48',
        'K33': '33',
        'K29': '29',
        'K27': '27',
        'K11': '11',
        'K6':  '6'
    }
    
    def number_to_letter(num):
        """Convert number to letter: 1->A, 2->B, 3->C, etc."""
        return chr(ord('A') + num - 1)
    
    return '-'.join(f"{number_to_letter(s)}{lysine_map.get(k, '?')}{number_to_letter(d)}" for s, k, d in edges)


# ==================================
# ==== Jeff all lysines Nomenclature Short Hand Conversion
# ==================================

def conjugated_lysines_to_chemical_all_node_nomenclature(conjugated_lysines):
    """
    Convert conjugated lysines format to Jeff's shorthand tree nomenclature system
    
    Args:
        conjugated_lysines (list): List of connections in format [[src, linkage, dst], ...]
                                  Example: [[1, 'K63', 2], [2, 'K48', 3], [1, 'K33', 4], [4, 'K11', 5]]
                                  Order matters - represents preorder traversal
    
    Returns:
        str: Jeff's shorthand nomenclature string
    
    Level System:
        Level 1 = A, Level 2 = B, Level 3 = C, Level 4 = D, etc.
    
    Position Mapping:
        K63: evens with uppercase letter (e.g., B2, C2, D4)
        K48: odds with uppercase letter (e.g., B1, C1, D3)  
        K33: evens with uppercase letter + ' (e.g., B'2, C'2, D'4)
        K29: odds with uppercase letter + ' (e.g., B'1, C'1, D'3)
        K11: evens with lowercase letter (e.g., b2, c2, d4)
        K6:  odds with lowercase letter (e.g., b1, c1, d3)
        K27: evens with lowercase letter + ' (e.g., b'2, c'2, d'4)
        M1:  odds with lowercase letter + ' (e.g., b'1, c'1, d'3)
    
    Formula: position = ((parent_letter_size + parent_number) - 1) * 2 + child_number
    
    Where:
        parent_letter_size: Numeric value based on parent's notation type
            - 0 for uppercase without prime (e.g., A1, B2)
            - 2 for uppercase with prime (e.g., A'1, B'2)  
            - 4 for lowercase without prime (e.g., a1, b2)
            - 6 for lowercase with prime (e.g., a'1, b'2)
        parent_number: The numeric part of the parent's notation (e.g., 2 from "B2" or "b'2")
        - 1 to count the number of nodes before you reach your node otherwise you also include your node
        child_number: Position increment based on lysine type
            - +1 for K48, K29, K6, M1 (odd positions)
            - +2 for K63, K33, K11, K27 (even positions)

        lysine_to_notation = {
        'K63': ('upper', False, 2),   # uppercase, no prime, child_number=+2
        'K48': ('upper', False, 1),   # uppercase, no prime, child_number=+1
        'K33': ('upper', True, 2),    # uppercase, prime, child_number=+2
        'K29': ('upper', True, 1),    # uppercase, prime, child_number=+1
        'K11': ('lower', False, 2),   # lowercase, no prime, child_number=+2
        'K6':  ('lower', False, 1),   # lowercase, no prime, child_number=+1
        'K27': ('lower', True, 2),    # lowercase, prime, child_number=+2
        'M1':  ('lower', True, 1)     # lowercase, prime, child_number=+1
    }
    """
    if not conjugated_lysines:
        return "A1"  # Single ubiquitin
    
    # Lysine to notation mapping with child_number increments
    # K63, K48, K33, K29, K27, K11, K6, M1 = 1, 2, 1, 2, 1, 2, 1, 2 
    lysine_to_notation = {
        'K63': ('upper', False, 2),   # uppercase, no prime, child_number=+2
        'K48': ('upper', False, 1),   # uppercase, no prime, child_number=+1
        'K33': ('upper', True, 2),    # uppercase, prime, child_number=+2
        'K29': ('upper', True, 1),    # uppercase, prime, child_number=+1
        'K11': ('lower', False, 2),   # lowercase, no prime, child_number=+2
        'K6':  ('lower', False, 1),   # lowercase, no prime, child_number=+1
        'K27': ('lower', True, 2),    # lowercase, prime, child_number=+2
        'M1':  ('lower', True, 1)     # lowercase, prime, child_number=+1
    }

    # Alternative mapping examples exactly - kept for reference
    lysine_to_notation_better = {
        'K63': ('upper', False, 1),   # uppercase, no prime, child_number=+2
        'K48': ('upper', False, 2),   # uppercase, no prime, child_number=+1
        'K33': ('upper', True, 1),    # uppercase, prime, child_number=+2
        'K29': ('upper', True, 2),    # uppercase, prime, child_number=+1
        'K27': ('lower', False, 1),    # lowercase, prime, child_number=+2
        'K11': ('lower', False, 2),   # lowercase, no prime, child_number=+2
        'K6':  ('lower', True, 1),   # lowercase, no prime, child_number=+1
        'M1':  ('lower', True, 2)     # lowercase, prime, child_number=+1
    }
    
    # Parent letter size mapping
    letter_size_mapping = {
        'upper_no_prime': 0,    # uppercase without prime
        'upper_prime': 2,       # uppercase with prime
        'lower_no_prime': 4,    # lowercase without prime
        'lower_prime': 6        # lowercase with prime
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
    
    # Assign shorthand positions
    nomenclature_map = {}  # node -> shorthand_notation
    
    def assign_shorthand_positions(node, level):
        if node == root:
            # Root node is always A1
            nomenclature_map[node] = "A1"
        else:
            # Get parent info and linkage to determine notation
            parent, linkage = parents[node]
            
            # Get notation rules for this linkage
            notation_info = lysine_to_notation.get(linkage)
            if notation_info is None:
                # Unknown linkage, skip
                return
            
            case_type, has_asterisk, child_number = notation_info
            
            # Determine level letter
            level_letter = chr(ord('A') + level)
            if case_type == 'lower':
                level_letter = level_letter.lower()
            
            # Calculate position based on formula: ((parent_letter_size + parent_number) - 1) * 2 + child_number
            if parent == root:
                # Direct children of root use the child_number based on lysine type
                position = child_number
            else:
                # Get parent's notation to extract parent_letter_size and parent_number
                parent_notation = nomenclature_map.get(parent, "")
                
                # Extract parent_number from notation
                import re
                match = re.search(r"(\d+)$", parent_notation)
                parent_number = int(match.group(1)) if match else 1
                
                # Determine parent_letter_size based on parent's notation type
                parent_letter_size = 0  # default
                if parent_notation:
                    if parent_notation[0].isupper():
                        if "'" in parent_notation:
                            parent_letter_size = 2  # uppercase with prime
                        else:
                            parent_letter_size = 0  # uppercase without prime
                    else:  # lowercase
                        if "'" in parent_notation:
                            parent_letter_size = 6  # lowercase with prime
                        else:
                            parent_letter_size = 4  # lowercase without prime
                
                # Apply the formula
                position = ((parent_letter_size + parent_number) - 1) * 2 + child_number
            
            # Build the notation string
            prime = "'" if has_asterisk else ""
            notation = f"{level_letter}{prime}{position}"
            nomenclature_map[node] = notation
        
        # Recursively assign positions to children (maintaining preorder)
        if node in children:
            # Sort children by their appearance order in the original list
            sorted_children = sorted(children[node], key=lambda x: x[2])  # Sort by order
            
            for child, linkage, _ in sorted_children:
                assign_shorthand_positions(child, level + 1)
    
    # Start assignment from root
    assign_shorthand_positions(root, 0)
    
    # Build nomenclature list by collecting all notations and reorder by level then position
    all_notations = []
    for node in all_nodes:
        if node in nomenclature_map:
            notation = nomenclature_map[node]
            
            # Parse notation to extract level letter and position number for sorting
            import re
            match = re.match(r"^([a-zA-Z])('?)(\d+)$", notation)
            if match:
                level_letter, prime, position_str = match.groups()
                level_num = ord(level_letter.upper()) - ord('A')
                position_num = int(position_str)
                
                all_notations.append((level_num, position_num, notation))
    
    # Build the final nomenclature string with comma separators
    nomenclature_parts = [notation for _, _, notation in all_notations]
    
    return ','.join(nomenclature_parts)


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



# ==================================
# ===== Deprecated Functions =======
# ==================================

# ==================================
# ==== K48-K63 Nomenclature Conversion (old)
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
# ==== all lysines Nomenclature Conversion
# ==================================

def conjugated_lysines_all_lysines_nomenclature(conjugated_lysines):
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