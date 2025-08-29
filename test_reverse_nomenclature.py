#!/usr/bin/env python3

import sys
from pathlib import Path

# Add the backend source to Python path
current_file = Path(__file__).resolve()
project_root = current_file.parent
sys.path.insert(0, str(project_root))
backend_src = project_root / 'back_end' / 'src'
sys.path.insert(0, str(backend_src))

from nomenclature import conjugated_lysines_to_nomenclature

def test_reverse_nomenclature():
    """Test the reverse conversion from conjugated lysines to nomenclature"""
    
    # Test case 1: Simple linear chain
    test_case_1 = [[1, 'K63', 2], [2, 'K63', 3]]
    expected_1 = "A1B1C1"
    result_1 = conjugated_lysines_to_nomenclature(test_case_1)
    print(f"Test 1: {test_case_1}")
    print(f"Expected: {expected_1}")
    print(f"Got: {result_1}")
    print(f"Pass: {result_1 == expected_1}\n")
    
    # Test case 2: Branching structure (provided example)
    test_case_2 = [[1, 'K63', 2], [2, 'K63', 3], [1, 'K48', 4], [4, 'K48', 5]]
    result_2 = conjugated_lysines_to_nomenclature(test_case_2)
    print(f"Test 2: {test_case_2}")
    print(f"Got: {result_2}")
    
    # Let's trace through this manually:
    # Node 1 is root -> A1
    # Node 2 connects via K63 from 1 -> position 1 (2*1-1=1) -> B1  
    # Node 4 connects via K48 from 1 -> position 2 (2*1=2) -> B2
    # Node 3 connects via K63 from 2 -> position 1 (2*1-1=1) -> C1
    # Node 5 connects via K48 from 4 -> position 4 (2*2=4) -> C4
    expected_2 = "A1B1B2C1C4"
    print(f"Expected: {expected_2}")
    print(f"Pass: {result_2 == expected_2}\n")
    
    # Test case 3: Single ubiquitin
    test_case_3 = []
    expected_3 = "A1"
    result_3 = conjugated_lysines_to_nomenclature(test_case_3)
    print(f"Test 3: {test_case_3}")
    print(f"Expected: {expected_3}")
    print(f"Got: {result_3}")
    print(f"Pass: {result_3 == expected_3}\n")
    
    # Test case 4: More complex branching with correct preorder traversal
    test_case_4 = [[1, 'K63', 2], [1, 'K48', 4], [2, 'K63', 3], [4, 'K48', 5]]
    result_4 = conjugated_lysines_to_nomenclature(test_case_4)
    print(f"Test 4: {test_case_4}")
    print(f"Got: {result_4}")
    # Expected: A1B1B2C1C4
    # 1 -> A1 (root)
    # 2 connects K63 from 1 -> B1 (first child via K63)
    # 4 connects K48 from 1 -> B2 (second child via K48)
    # 3 connects K63 from 2 -> C1 (child of B1)
    # 5 connects K48 from 4 -> C4 (child of B2)
    expected_4 = "A1B1B2C1C4"
    print(f"Expected: {expected_4}")
    print(f"Pass: {result_4 == expected_4}\n")

    # Test case 4: More complex branching with correct preorder traversal
    test_case_4 = [[1, 'K63', 2], [2, 'K63', 3], [1, 'K48', 4], [4, 'K48', 5]]
    result_4 = conjugated_lysines_to_nomenclature(test_case_4)
    print(f"Test 4: {test_case_4}")
    print(f"Got: {result_4}")
    # Expected: A1B1B2C1C4
    # 1 -> A1 (root)
    # 2 connects K63 from 1 -> B1 (first child via K63)
    # 4 connects K48 from 1 -> B2 (second child via K48)
    # 3 connects K63 from 2 -> C1 (child of B1)
    # 5 connects K48 from 4 -> C4 (child of B2)
    expected_4 = "A1B1B2C1C4"
    print(f"Expected: {expected_4}")
    print(f"Pass: {result_4 == expected_4}\n")

if __name__ == "__main__":
    test_reverse_nomenclature()
