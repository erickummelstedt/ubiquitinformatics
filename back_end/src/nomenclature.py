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



