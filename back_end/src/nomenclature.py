from pathlib import Path
import pandas as pd
import sys

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
Level System:

A = level 1, B = level 2, C = level 3, D = level 4, etc.
Numbering System:

Odd numbers (1, 3, 5, 7, ...) = K63 linkages
Even numbers (2, 4, 6, 8, ...) = K48 linkages
Level Capacity:

Level A: n = 1 (so A1 only)
Level B: n = 2 (so B1, B2)
Level C: n = 4 (so C1, C2, C3, C4)
Level D: n = 8 (so D1, D2, D3, D4, D5, D6, D7, D8)
Pattern: each level has 2^(level-1) positions
Binding Rules:

Each position can only bind to specific positions in the next level
The rule is: position x can bind to positions (2x-1) and (2x) in the next level
Examples:
A1 → B1, B2
B1 → C1, C2
B2 → C3, C4
C4 → D7, D8
"""


def parse_tree_nomenclature(tree_string):
    """
    Parse tree nomenclature string into individual positions
    Example: 'A1B2C4D6' -> [('A', 1), ('B', 2), ('C', 4), ('D', 6)]
    """
    positions = []
    i = 0
    while i < len(tree_string):
        if tree_string[i].isalpha():
            level = tree_string[i]
            i += 1
            # Extract the number
            num_str = ''
            while i < len(tree_string) and tree_string[i].isdigit():
                num_str += tree_string[i]
                i += 1
            if num_str:
                positions.append((level, int(num_str)))
        else:
            i += 1
    return positions

def validate_tree_structure(positions):
    """
    Validate that the tree structure follows the binding rules
    Position x can only bind to positions (2x-1) and (2x) in the next level
    """
    for i in range(len(positions) - 1):
        current_level, current_pos = positions[i]
        next_level, next_pos = positions[i + 1]
        
        # Check if the binding rule is satisfied
        expected_pos_1 = 2 * current_pos - 1
        expected_pos_2 = 2 * current_pos
        
        if next_pos not in [expected_pos_1, expected_pos_2]:
            return False, f"Invalid binding: {current_level}{current_pos} cannot bind to {next_level}{next_pos}"
    
    return True, "Valid tree structure"

def determine_linkage_type(position):
    """
    Determine linkage type based on position number
    Odd = K63, Even = K48
    """
    return "K63" if position % 2 == 1 else "K48"

def tree_to_graphical(tree_string):
    """
    Convert tree nomenclature to graphical description with preorder numbering
    """
    positions = parse_tree_nomenclature(tree_string)
    
    # Validate structure
    is_valid, message = validate_tree_structure(positions)
    if not is_valid:
        return None, message
    
    # Convert to graphical representation
    graphical_chain = []
    
    for i in range(len(positions) - 1):
        current_node = i + 1  # Preorder numbering starts at 1
        next_node = i + 2
        
        # Linkage type determined by target position
        _, target_pos = positions[i + 1]
        linkage_type = determine_linkage_type(target_pos)
        
        graphical_chain.append(f"{current_node} → {linkage_type} → {next_node}")
    
    return graphical_chain, "Success"

def convert_ubiquitin_nomenclature(tree_string):
    """
    Single function to convert tree nomenclature to graphical representation
    Handles both linear chains and branching structures with proper preorder numbering
    
    Args:
        tree_string (str): Tree nomenclature string (e.g., 'A1B2C4D6' or 'A1B1B2C3D5')
    
    Returns:
        dict: {
            'success': bool,
            'result': list or str,  # graphical chain if success, error message if not
            'input': str            # original input for reference
        }
    """
    try:
        # Parse the input
        positions = []
        i = 0
        while i < len(tree_string):
            if tree_string[i].isalpha():
                level = tree_string[i]
                i += 1
                # Extract the number
                num_str = ''
                while i < len(tree_string) and tree_string[i].isdigit():
                    num_str += tree_string[i]
                    i += 1
                if num_str:
                    positions.append((level, int(num_str)))
            else:
                i += 1
        
        if not positions:
            return {
                'success': False,
                'result': 'Invalid input: No valid positions found',
                'input': tree_string
            }
        
        # Group positions by level
        levels = {}
        for level, pos in positions:
            if level not in levels:
                levels[level] = []
            levels[level].append(pos)
        
        # Sort positions within each level
        for level in levels:
            levels[level].sort()
        
        # Build tree structure
        tree = {}
        level_names = sorted(levels.keys())
        
        # Create parent-child relationships
        for i in range(len(level_names) - 1):
            current_level = level_names[i]
            next_level = level_names[i + 1]
            
            for current_pos in levels[current_level]:
                expected_pos_1 = 2 * current_pos - 1
                expected_pos_2 = 2 * current_pos
                
                children = []
                for next_pos in levels[next_level]:
                    if next_pos in [expected_pos_1, expected_pos_2]:
                        children.append(next_pos)
                
                if children:
                    tree[(current_level, current_pos)] = [(next_level, child) for child in children]
        
        # Validate that all positions have valid parents (except root)
        for i in range(1, len(positions)):
            current_level, current_pos = positions[i]
            
            # Find potential parents
            parent_found = False
            for j in range(i):
                parent_level, parent_pos = positions[j]
                expected_pos_1 = 2 * parent_pos - 1
                expected_pos_2 = 2 * parent_pos
                
                if current_pos in [expected_pos_1, expected_pos_2]:
                    if ord(parent_level) == ord(current_level) - 1:
                        parent_found = True
                        break
            
            if not parent_found:
                return {
                    'success': False,
                    'result': f"Validation: Invalid binding: {current_level}{current_pos} has no valid parent",
                    'input': tree_string
                }
        
        # Assign preorder numbers using DFS traversal
        def assign_preorder_numbers(node, node_counter):
            position_to_node[node] = node_counter[0]
            node_counter[0] += 1
            
            # Visit children in order
            if node in tree:
                for child in sorted(tree[node], key=lambda x: x[1]):  # Sort by position number
                    assign_preorder_numbers(child, node_counter)
        
        position_to_node = {}
        node_counter = [1]  # Use list to allow modification in nested function
        
        # Start DFS from root
        root = positions[0]
        assign_preorder_numbers(root, node_counter)
        
        # Build graphical chain
        graphical_chain = []
        
        def build_connections(node):
            if node in tree:
                for child in sorted(tree[node], key=lambda x: x[1]):
                    parent_node_num = position_to_node[node]
                    child_node_num = position_to_node[child]
                    
                    # Linkage type determined by target position
                    _, child_pos = child
                    linkage_type = "K63" if child_pos % 2 == 1 else "K48"
                    
                    graphical_chain.append(f"{parent_node_num} → {linkage_type} → {child_node_num}")
                    
                    # Recursively build connections for children
                    build_connections(child)
        
        build_connections(root)
        
        return {
            'success': True,
            'result': graphical_chain,
            'input': tree_string
        }
        
    except Exception as e:
        return {
            'success': False,
            'result': f"Error processing input: {str(e)}",
            'input': tree_string
        }


def build_polyubiquitin_from_nomenclature(tree_string):
    """
    Convert tree nomenclature to actual polyubiquitin structure using ubiquitin_building_all
    
    Args:
        tree_string (str): Tree nomenclature string (e.g., 'A1B1B2C3D5')
    
    Returns:
        dict: The final polyubiquitin structure or error information
    """
    # First convert the nomenclature to graphical representation
    result = convert_ubiquitin_nomenclature(tree_string)
    
    connections = result['result']
    
    # Start with the base ubiquitin molecule
    current_structure = histag_ubi_ubq_1
    
    # Apply each connection iteratively
    for connection in connections:
        # Parse the connection string: "1 → K63 → 2"
        parts = connection.split(' → ')
        ubiquitin_number = int(parts[0])
        lysine_residue = parts[1]  # K63 or K48
        # parts[2] is the target ubiquitin number (not needed for ubiquitin_building_all)
        
        current_structure, context = ubiquitin_building_all(
            parent_dictionary=current_structure,
            ubi_molecule_to_add=ubi_ubq_1,
            ubiquitin_number=ubiquitin_number,
            lysine_residue=lysine_residue
        )

    output_structure, output_context = iterate_through_ubiquitin(current_structure)
    
    return output_structure, connections
