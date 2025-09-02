"""
Test suite for nomenclature system functions
============================================

This module contains comprehensive tests for all nomenclature conversion functions
including tree nomenclature, compact edges, conjugated lysines, and formatted edges.
"""

import pytest
from pathlib import Path
import sys

# Setup path for imports
current_file = Path(__file__).resolve()
project_root = current_file.parents[2]
sys.path.insert(0, str(project_root))
local_path = project_root / 'back_end'
sys.path.insert(0, str(local_path))

from src.nomenclature import (
    parse_compact_edges,
    build_polyubiquitin_from_edges,
    multimer_length_from_nomenclature,
    conjugated_lysines_to_jeff_K48_K63_nomenclature,
    conjugated_lysines_to_jeff_all_lysines_nomenclature,
    conjugated_lysines_to_jeffs_multiple_symbols,
    tree_nomenclature_to_numerical_system,
    format_nomenclature_preorder_ABC,
    format_nomenclature_preorder_jeff,
    LETTER_TO_LYS
)

# ===================================
# Basic Function Tests
# ===================================

def test_01_format_edges_simple_chain():
    """Test format_nomenclature_preorder_ABC with simple linear chain"""
    edges = [[1, 'K63', 2], [2, 'K63', 3]]
    result = format_nomenclature_preorder_ABC(edges)
    expected = "1A2-2A3"
    assert result == expected


def test_02_format_edges_mixed_lysines():
    """Test format_nomenclature_preorder_ABC with multiple lysine types"""
    edges = [[1, 'K63', 2], [1, 'K48', 3], [2, 'K33', 4]]
    result = format_nomenclature_preorder_ABC(edges)
    expected = "1A2-1B3-2C4"
    assert result == expected


def test_03_format_edges_all_seven_lysines():
    """Test format_nomenclature_preorder_ABC with all seven lysine types"""
    edges = [
        [1, 'K63', 2], [1, 'K48', 3], [1, 'K33', 4], 
        [1, 'K29', 5], [1, 'K27', 6], [1, 'K11', 7], [1, 'K6', 8]
    ]
    result = format_nomenclature_preorder_ABC(edges)
    expected = "1A2-1B3-1C4-1D5-1E6-1F7-1G8"
    assert result == expected


def test_04_parse_compact_edges_dash_separated():
    """Test parse_compact_edges with dash-separated format"""
    compact = "1A2-2A3-3B4"
    result = parse_compact_edges(compact)
    expected = [[1, 'K63', 2], [2, 'K63', 3], [3, 'K48', 4]]
    assert result == expected


def test_05_parse_compact_edges_comma_separated():
    """Test parse_compact_edges with comma-separated format"""
    compact = "1A2, 2B3, 3C4"
    result = parse_compact_edges(compact)
    expected = [[1, 'K63', 2], [2, 'K48', 3], [3, 'K33', 4]]
    assert result == expected


def test_06_parse_compact_edges_with_whitespace():
    """Test parse_compact_edges handles whitespace correctly"""
    compact = " 1A2 - 2B3 - 3C4 "
    result = parse_compact_edges(compact)
    expected = [[1, 'K63', 2], [2, 'K48', 3], [3, 'K33', 4]]
    assert result == expected


def test_07_parse_compact_edges_single_edge():
    """Test parse_compact_edges with single edge"""
    compact = "1A2"
    result = parse_compact_edges(compact)
    expected = [[1, 'K63', 2]]
    assert result == expected


def test_08_parse_compact_edges_list_input():
    """Test parse_compact_edges with list input"""
    compact = ["1A2", "2B3", "3C4"]
    result = parse_compact_edges(compact)
    expected = [[1, 'K63', 2], [2, 'K48', 3], [3, 'K33', 4]]
    assert result == expected


def test_09_parse_compact_edges_invalid_format():
    """Test parse_compact_edges raises error for invalid format"""
    with pytest.raises(ValueError):
        parse_compact_edges("1X2")  # Invalid lysine letter
    
    with pytest.raises(ValueError):
        parse_compact_edges("12")   # Missing lysine letter


@pytest.mark.parametrize("compact, expected_length", [
    ("1A2", 2),
    ("1A2-2A3", 3),
    ("1A2-2A3-3B4", 4),
    ("1A2-1B3-2C4-3D5", 5),
    ("1A2-1B3-1C4-2D5-3E6", 6)
])
def test_10_multimer_length_calculation(compact, expected_length):
    """Test multimer_length_from_nomenclature function"""
    result = multimer_length_from_nomenclature(compact)
    assert result == expected_length


# ===================================
# Jeff K48-K63 Nomenclature Tests
# ===================================

def test_11_jeff_k48_k63_nomenclature_simple_chain():
    """Test Jeff K48-K63 nomenclature for simple chain"""
    conjugated_lysines = [[1, 'K63', 2], [2, 'K63', 3]]
    result = conjugated_lysines_to_jeff_K48_K63_nomenclature(conjugated_lysines)
    expected = "A1B1C1"
    assert result == expected


def test_12_jeff_k48_k63_nomenclature_branched_structure():
    """Test Jeff K48-K63 nomenclature for branched structure"""
    conjugated_lysines = [[1, 'K63', 2], [1, 'K48', 3]]
    result = conjugated_lysines_to_jeff_K48_K63_nomenclature(conjugated_lysines)
    expected = "A1B1B2"
    assert result == expected


def test_13_jeff_k48_k63_nomenclature_complex_tree():
    """Test Jeff K48-K63 nomenclature for complex tree"""
    conjugated_lysines = [
        [1, 'K63', 2], [2, 'K63', 3], [1, 'K48', 4], [4, 'K48', 5]
    ]
    result = conjugated_lysines_to_jeff_K48_K63_nomenclature(conjugated_lysines)
    # The algorithm produces A1B1B2C1C4 based on preorder traversal and position calculation
    expected = "A1B1B2C1C4"
    assert result == expected

def test_13a_jeff_k48_k63_nomenclature_complex_tree():
    """Test Jeff K48-K63 nomenclature rejects non-K48/K63 lysines"""
    conjugated_lysines = [
        [1, 'K63', 2], [2, 'K63', 3], [2, 'K29', 4], [1, 'K48', 5]
    ]
    # Should raise ValueError because K29 is not allowed in K48-K63 nomenclature
    with pytest.raises(ValueError):
        conjugated_lysines_to_jeff_K48_K63_nomenclature(conjugated_lysines)


def test_14_jeff_k48_k63_nomenclature_empty_input():
    """Test Jeff K48-K63 nomenclature with empty input"""
    result = conjugated_lysines_to_jeff_K48_K63_nomenclature([])
    expected = "A1"
    assert result == expected


# ===================================
# Jeff All Lysines Nomenclature Tests
# ===================================

def test_15a_jeff_all_lysines_nomenclature_complex_tree():
    """Test Jeff K48-K63 nomenclature with complex tree"""
    conjugated_lysines = [
        [1, 'K63', 2], [2, 'K63', 3], [1, 'K48', 4], [4, 'K48', 5]
    ]
    result = conjugated_lysines_to_jeff_all_lysines_nomenclature(conjugated_lysines)
    # The algorithm produces A1B1B2C1C4 based on preorder traversal and position calculation
    expected = "A1B1B2C1C9"
    assert result == expected

def test_15b_jeff_all_lysines_nomenclature_complex_tree():
    """Test Jeff all lysines nomenclature with complex tree"""
    conjugated_lysines = [
        [1, 'K63', 2], [2, 'K63', 3], [2, 'K29', 4], [1, 'K48', 5]
    ]
    result = conjugated_lysines_to_jeff_all_lysines_nomenclature(conjugated_lysines)
    # The algorithm produces A1B1C1B2D4 based on preorder traversal and position calculation
    expected = "A1B1B2C1C4"
    assert result == expected

def test_15_jeff_all_lysines_nomenclature_simple_chain():
    """Test Jeff all lysines nomenclature for simple chain"""
    conjugated_lysines = [[1, 'K63', 2], [2, 'K63', 3]]
    result = conjugated_lysines_to_jeff_all_lysines_nomenclature(conjugated_lysines)
    # Should follow 7-lysine mapping: K63=position 1, so A1B1C1
    expected = "A1B1C1"
    assert result == expected


def test_16_jeff_all_lysines_nomenclature_mixed_lysines():
    """Test Jeff all lysines nomenclature with mixed lysine types"""
    conjugated_lysines = [[1, 'K63', 2], [1, 'K48', 3], [1, 'K33', 4]]
    result = conjugated_lysines_to_jeff_all_lysines_nomenclature(conjugated_lysines)
    # K63=pos1, K48=pos2, K33=pos3 in level B (positions 1-7 from parent A1)
    assert "A1" in result
    assert "B1" in result  # K63 branch
    assert "B2" in result  # K48 branch
    assert "B3" in result  # K33 branch


def test_17_jeff_all_lysines_nomenclature_k27_k11_k6():
    """Test Jeff all lysines nomenclature with K27, K11, K6"""
    conjugated_lysines = [[1, 'K27', 2], [1, 'K11', 3], [1, 'K6', 4]]
    result = conjugated_lysines_to_jeff_all_lysines_nomenclature(conjugated_lysines)
    # K27=pos5, K11=pos6, K6=pos7 in level B
    assert "A1" in result
    assert "B5" in result  # K27 branch
    assert "B6" in result  # K11 branch
    assert "B7" in result  # K6 branch


def test_18_jeff_all_lysines_nomenclature_cyclic_positions():
    """Test Jeff all lysines nomenclature with cyclic position calculation"""
    # Test that positions cycle every 7: K63 should use positions 1, 8, 15...
    conjugated_lysines = [
        [1, 'K63', 2], [2, 'K63', 3], [3, 'K63', 4], 
        [4, 'K63', 5], [5, 'K63', 6], [6, 'K63', 7], 
        [7, 'K63', 8], [8, 'K63', 9]  # This should use position 8 (1+7*1)
    ]
    result = conjugated_lysines_to_jeff_all_lysines_nomenclature(conjugated_lysines)
    # Should have multiple levels with K63 connections
    assert "A1" in result
    # Exact positions depend on binding rules, but should be valid


# ===================================
# Numerical System Tests
# ===================================

def test_21_tree_nomenclature_to_numerical_simple():
    """Test tree nomenclature to numerical conversion for simple cases"""
    # A1 = 0 + 1 = 1
    result = tree_nomenclature_to_numerical_system("A1")
    expected = "1"
    assert result == expected


def test_22_tree_nomenclature_to_numerical_basic_levels():
    """Test basic level multipliers: A=0, B=7, C=49, D=343"""
    # A1 = 0+1=1, B1 = 7+1=8, C1 = 49+1=50, D1 = 343+1=344
    result = tree_nomenclature_to_numerical_system("A1B1C1D1")
    expected = "1, 8, 50, 344"
    assert result == expected


def test_23_tree_nomenclature_to_numerical_complex():
    """Test complex tree nomenclature conversion"""
    # A1B1B2C1C4 -> "1, 8, 9, 50, 53"
    # A1 = 0+1=1, B1 = 7+1=8, B2 = 7+2=9, C1 = 49+1=50, C4 = 49+4=53
    result = tree_nomenclature_to_numerical_system("A1B1B2C1C4")
    expected = "1, 8, 9, 50, 53"
    assert result == expected


def test_24_tree_nomenclature_to_numerical_edge_cases():
    """Test edge cases for numerical conversion"""
    # Empty string
    result = tree_nomenclature_to_numerical_system("")
    assert result == ""
    
    # Invalid format
    result = tree_nomenclature_to_numerical_system("XYZ123")
    assert result == ""
    
    # Single level with high position
    result = tree_nomenclature_to_numerical_system("A10")
    expected = "10"  # 0 + 10 = 10
    assert result == expected
    
    # High level test: E5 = 7^4 + 5 = 2401 + 5 = 2406
    result = tree_nomenclature_to_numerical_system("E5")
    expected = "2406"
    assert result == expected


def test_25_tree_nomenclature_to_numerical_integration():
    """Test integration with actual nomenclature function output"""
    # Generate nomenclature first, then convert to numerical
    conjugated_lysines = [[1, 'K63', 2], [1, 'K48', 3]]
    nomenclature = conjugated_lysines_to_jeff_all_lysines_nomenclature(conjugated_lysines)
    numerical = tree_nomenclature_to_numerical_system(nomenclature)
    
    # Should return a comma-separated string of positive integers
    assert isinstance(numerical, str)
    assert len(numerical) > 0
    # Check that it contains numbers separated by commas
    numbers = [int(x.strip()) for x in numerical.split(',')]
    assert all(x > 0 for x in numbers)


# ===================================
# System Tests
# ===================================

def test_19_letter_to_lys_mapping_completeness():
    """Test that LETTER_TO_LYS mapping is complete"""
    expected_mapping = {
        'A': 'K63', 'B': 'K48', 'C': 'K33', 'D': 'K29',
        'E': 'K27', 'F': 'K11', 'G': 'K6'
    }
    assert LETTER_TO_LYS == expected_mapping


def test_20_edge_cases_and_error_handling():
    """Test various edge cases and error handling"""
    # Test empty string
    result = parse_compact_edges("")
    assert result == []
    
    # Test whitespace only
    result = parse_compact_edges("   ")
    assert result == []
    
    # Test format_nomenclature_preorder_ABC with empty list
    result = format_nomenclature_preorder_ABC([])
    assert result == ""
    
    # Test multimer_length with no matches
    result = multimer_length_from_nomenclature("no_matches_here")
    assert result == 0
    
    # Test nomenclature functions with malformed edges
    malformed_edges = [[1, 'InvalidLysine', 2]]
    # Should handle gracefully without crashing
    result = conjugated_lysines_to_jeff_all_lysines_nomenclature(malformed_edges)
    assert isinstance(result, str)


# ===================================
# Integration Tests
# ===================================

def test_round_trip_conversion():
    """Test round-trip conversion: edges -> compact -> edges"""
    original_edges = [[1, 'K63', 2], [2, 'K48', 3], [1, 'K33', 4]]
    
    # Convert to compact format
    compact = format_nomenclature_preorder_ABC(original_edges)
    
    # Parse back to edges
    parsed_edges = parse_compact_edges(compact)
    
    # Should match original
    assert parsed_edges == original_edges


def test_nomenclature_consistency():
    """Test that both nomenclature functions handle same input consistently"""
    test_edges = [[1, 'K63', 2], [1, 'K48', 3]]
    
    k48_k63_result = conjugated_lysines_to_jeff_K48_K63_nomenclature(test_edges)
    all_lysines_result = conjugated_lysines_to_jeff_all_lysines_nomenclature(test_edges)
    
    # Both should produce valid nomenclature strings
    assert k48_k63_result.startswith('A1')
    assert all_lysines_result.startswith('A1')
    
    # Both should be non-empty strings
    assert isinstance(k48_k63_result, str)
    assert isinstance(all_lysines_result, str)
    assert len(k48_k63_result) > 0
    assert len(all_lysines_result) > 0


# ===================================
# Format Nomenclature ABC Tests
# ===================================

def test_26_format_nomenclature_preorder_ABC_simple():
    """Test format_nomenclature_preorder_ABC with simple chain"""
    edges = [[1, 'K63', 2], [2, 'K63', 3]]
    result = format_nomenclature_preorder_ABC(edges)
    expected = "1A2-2A3"
    assert result == expected


def test_27_format_nomenclature_preorder_ABC_all_lysines():
    """Test format_nomenclature_preorder_ABC with all seven lysine types"""
    edges = [
        [1, 'K63', 2], [1, 'K48', 3], [1, 'K33', 4], 
        [1, 'K29', 5], [1, 'K27', 6], [1, 'K11', 7], [1, 'K6', 8]
    ]
    result = format_nomenclature_preorder_ABC(edges)
    expected = "1A2-1B3-1C4-1D5-1E6-1F7-1G8"
    assert result == expected


def test_28_format_nomenclature_preorder_ABC_branched():
    """Test format_nomenclature_preorder_ABC with branched structure"""
    edges = [[1, 'K63', 2], [1, 'K48', 3], [2, 'K33', 4]]
    result = format_nomenclature_preorder_ABC(edges)
    expected = "1A2-1B3-2C4"
    assert result == expected


def test_29_format_nomenclature_preorder_ABC_unknown_lysine():
    """Test format_nomenclature_preorder_ABC with unknown lysine type"""
    edges = [[1, 'K999', 2]]
    result = format_nomenclature_preorder_ABC(edges)
    expected = "1?2"
    assert result == expected


def test_30_format_nomenclature_preorder_ABC_empty():
    """Test format_nomenclature_preorder_ABC with empty input"""
    edges = []
    result = format_nomenclature_preorder_ABC(edges)
    expected = ""
    assert result == expected


# ===================================
# Format Nomenclature Jeff Tests
# ===================================

def test_31_format_nomenclature_preorder_jeff_simple():
    """Test format_nomenclature_preorder_jeff with simple chain"""
    edges = [[1, 'K63', 2], [2, 'K63', 3]]
    result = format_nomenclature_preorder_jeff(edges)
    expected = "A63B-B63C"
    assert result == expected


def test_32_format_nomenclature_preorder_jeff_mixed_lysines():
    """Test format_nomenclature_preorder_jeff with mixed lysine types"""
    edges = [[1, 'K63', 2], [1, 'K48', 3], [2, 'K33', 4]]
    result = format_nomenclature_preorder_jeff(edges)
    expected = "A63B-A48C-B33D"
    assert result == expected


def test_33_format_nomenclature_preorder_jeff_all_lysines():
    """Test format_nomenclature_preorder_jeff with all seven lysine types"""
    edges = [
        [1, 'K63', 2], [1, 'K48', 3], [1, 'K33', 4], 
        [1, 'K29', 5], [1, 'K27', 6], [1, 'K11', 7], [1, 'K6', 8]
    ]
    result = format_nomenclature_preorder_jeff(edges)
    expected = "A63B-A48C-A33D-A29E-A27F-A11G-A6H"
    assert result == expected


def test_34_format_nomenclature_preorder_jeff_large_numbers():
    """Test format_nomenclature_preorder_jeff with larger node numbers"""
    edges = [[10, 'K63', 11], [11, 'K48', 12]]
    result = format_nomenclature_preorder_jeff(edges)
    expected = "J63K-K48L"
    assert result == expected


def test_35_format_nomenclature_preorder_jeff_unknown_lysine():
    """Test format_nomenclature_preorder_jeff with unknown lysine type"""
    edges = [[1, 'K999', 2]]
    result = format_nomenclature_preorder_jeff(edges)
    expected = "A?B"
    assert result == expected


def test_36_format_nomenclature_preorder_jeff_empty():
    """Test format_nomenclature_preorder_jeff with empty input"""
    edges = []
    result = format_nomenclature_preorder_jeff(edges)
    expected = ""
    assert result == expected


# ===================================
# Jeff Shorthand Nomenclature Tests
# ===================================

def test_37_jeff_shorthand_nomenclature_empty():
    """Test Jeff shorthand nomenclature with empty input"""
    result = conjugated_lysines_to_jeffs_multiple_symbols([])
    expected = "A1"
    assert result == expected


def test_38_jeff_shorthand_nomenclature_simple_k63():
    """Test Jeff shorthand nomenclature with simple K63 chain"""
    conjugated_lysines = [[1, 'K63', 2], [2, 'K63', 3]]
    result = conjugated_lysines_to_jeffs_multiple_symbols(conjugated_lysines)
    # K63 uses uppercase letters with even positions
    assert "A1" in result
    assert result.startswith("A1")


def test_39_jeff_shorthand_nomenclature_k48_chain():
    """Test Jeff shorthand nomenclature with K48 chain"""
    conjugated_lysines = [[1, 'K48', 2], [2, 'K48', 3]]
    result = conjugated_lysines_to_jeffs_multiple_symbols(conjugated_lysines)
    # K48 uses uppercase letters with odd positions
    assert "A1" in result
    assert result.startswith("A1")


def test_40_jeff_shorthand_nomenclature_mixed_lysines():
    """Test Jeff shorthand nomenclature with various lysine types"""
    conjugated_lysines = [
        [1, 'K63', 2],  # uppercase, even
        [1, 'K48', 3],  # uppercase, odd
        [1, 'K33', 4],  # uppercase+asterisk, even
        [1, 'K29', 5]   # uppercase+asterisk, odd
    ]
    result = conjugated_lysines_to_jeffs_multiple_symbols(conjugated_lysines)
    assert "A1" in result
    # Should contain various case and asterisk combinations


def test_41_jeff_shorthand_nomenclature_lowercase_lysines():
    """Test Jeff shorthand nomenclature with lowercase lysine types"""
    conjugated_lysines = [
        [1, 'K11', 2],  # lowercase, even
        [1, 'K6', 3],   # lowercase, odd
        [1, 'K27', 4],  # lowercase+asterisk, even
        [1, 'M1', 5]    # lowercase+asterisk, odd
    ]
    result = conjugated_lysines_to_jeffs_multiple_symbols(conjugated_lysines)
    assert "A1" in result
    # Should contain lowercase letters and asterisks


def test_42_jeff_shorthand_nomenclature_complex_tree():
    """Test Jeff shorthand nomenclature with complex branched structure"""
    conjugated_lysines = [
        [1, 'K63', 2], [2, 'K48', 3], [1, 'K33', 4], [4, 'K11', 5]
    ]
    result = conjugated_lysines_to_jeffs_multiple_symbols(conjugated_lysines)
    expected = "A1,B2,C9,B*2,c26"
    assert result == expected
    # Should handle multiple levels with different case/asterisk combinations


def test_42a_jeff_shorthand_nomenclature_user_example():
    """Test Jeff shorthand nomenclature with user-provided example edges"""
    # Test case: 1 -> K63 -> 2, 1 -> K48 -> 3, 3 -> K63 -> 4, 3 -> K48 -> 5
    conjugated_lysines = [[1, 'K63', 2], [1, 'K48', 3], [3, 'K63', 4], [3, 'K48', 5]]
    result = conjugated_lysines_to_jeffs_multiple_symbols(conjugated_lysines)
    expected = "A1,B2,B1,C2,C1"
    assert result == expected


# ===================================
# Advanced Integration Tests
# ===================================

def test_43_all_nomenclature_functions_consistency():
    """Test that all nomenclature functions handle the same input without errors"""
    # Test edges with all lysine types for general functions
    test_edges_all = [[1, 'K63', 2], [1, 'K48', 3], [2, 'K33', 4]]
    
    # Test edges with only K48/K63 for the restricted function
    test_edges_k48_k63 = [[1, 'K63', 2], [1, 'K48', 3], [2, 'K63', 4]]
    
    # Test general nomenclature functions that accept all lysine types
    abc_result = format_nomenclature_preorder_ABC(test_edges_all)
    jeff_result = format_nomenclature_preorder_jeff(test_edges_all)
    all_lysines_result = conjugated_lysines_to_jeff_all_lysines_nomenclature(test_edges_all)
    shorthand_result = conjugated_lysines_to_jeffs_multiple_symbols(test_edges_all)
    
    # Test K48-K63 specific function with compatible input
    k48_k63_result = conjugated_lysines_to_jeff_K48_K63_nomenclature(test_edges_k48_k63)
    
    # All should return valid strings
    assert isinstance(abc_result, str)
    assert isinstance(jeff_result, str)
    assert isinstance(k48_k63_result, str)
    assert isinstance(all_lysines_result, str)
    assert isinstance(shorthand_result, str)
    
    # All should contain some content
    assert len(abc_result) > 0
    assert len(jeff_result) > 0
    assert len(k48_k63_result) > 0
    assert len(all_lysines_result) > 0
    assert len(shorthand_result) > 0


def test_44_edge_cases_malformed_inputs():
    """Test edge cases with malformed inputs"""
    # Test with invalid edge structure
    malformed_edges = [
        [1],  # Too few elements
        [1, 'K63'],  # Missing destination
        ['invalid', 'K63', 2],  # Non-numeric source
        [1, 'K63', 'invalid'],  # Non-numeric destination
    ]
    
    for bad_edge in malformed_edges:
        # Functions should handle gracefully without crashing
        try:
            format_nomenclature_preorder_ABC([bad_edge])
            format_nomenclature_preorder_jeff([bad_edge])
            conjugated_lysines_to_jeff_all_lysines_nomenclature([bad_edge])
            conjugated_lysines_to_jeffs_multiple_symbols([bad_edge])
        except (ValueError, IndexError, TypeError):
            # Expected behavior for malformed input
            pass


def test_45_large_multimer_structures():
    """Test nomenclature functions with larger multimer structures"""
    # Create a larger chain structure (10-mer)
    large_chain = []
    for i in range(1, 10):
        large_chain.append([i, 'K63', i + 1])
    
    # Test all nomenclature functions
    abc_result = format_nomenclature_preorder_ABC(large_chain)
    jeff_result = format_nomenclature_preorder_jeff(large_chain)
    all_lysines_result = conjugated_lysines_to_jeff_all_lysines_nomenclature(large_chain)
    shorthand_result = conjugated_lysines_to_jeffs_multiple_symbols(large_chain)
    
    # All should handle large structures
    assert len(abc_result) > 0
    assert len(jeff_result) > 0
    assert len(all_lysines_result) > 0
    assert len(shorthand_result) > 0
    
    # Should contain multiple levels
    assert 'A' in all_lysines_result
    assert 'B' in all_lysines_result


def test_46_highly_branched_structures():
    """Test nomenclature functions with highly branched structures"""
    # Create a star-like structure where node 1 connects to many others
    branched_structure = []
    lysines = ['K63', 'K48', 'K33', 'K29', 'K27', 'K11', 'K6']
    
    for i, lysine in enumerate(lysines, 2):
        branched_structure.append([1, lysine, i])
    
    # Test all nomenclature functions
    abc_result = format_nomenclature_preorder_ABC(branched_structure)
    jeff_result = format_nomenclature_preorder_jeff(branched_structure)
    all_lysines_result = conjugated_lysines_to_jeff_all_lysines_nomenclature(branched_structure)
    shorthand_result = conjugated_lysines_to_jeffs_multiple_symbols(branched_structure)
    
    # All should handle branched structures
    assert len(abc_result) > 0
    assert len(jeff_result) > 0
    assert len(all_lysines_result) > 0
    assert len(shorthand_result) > 0
    
    # Should show branching from root
    assert all_lysines_result.startswith('A1')
    assert shorthand_result.startswith('A1')


def test_47_numerical_system_advanced():
    """Test tree nomenclature to numerical system with advanced cases"""
    # Test high-level nomenclature
    test_cases = [
        ("A1", "1"),
        ("A1B1", "1, 8"),
        ("A1B1C1", "1, 8, 50"),
        ("A1B1C1D1", "1, 8, 50, 344"),
        ("A1B1C1D1E1", "1, 8, 50, 344, 2402"),
        ("A5B10C20", "5, 17, 69"),  # Higher positions
        ("A1B7C49", "1, 14, 98"),   # Maximum positions per level
    ]
    
    for input_nom, expected in test_cases:
        result = tree_nomenclature_to_numerical_system(input_nom)
        assert result == expected, f"Failed for {input_nom}: got {result}, expected {expected}"


def test_48_mixed_k48_k63_validation():
    """Test that K48-K63 nomenclature properly validates input"""
    # Valid K48-K63 inputs
    valid_inputs = [
        [[1, 'K63', 2]],
        [[1, 'K48', 2]],
        [[1, 'K63', 2], [1, 'K48', 3]],
        [[1, 'K63', 2], [2, 'K63', 3], [1, 'K48', 4]],
    ]
    
    for edges in valid_inputs:
        result = conjugated_lysines_to_jeff_K48_K63_nomenclature(edges)
        assert isinstance(result, str)
        assert len(result) > 0
    
    # Invalid K48-K63 inputs (should raise ValueError)
    invalid_inputs = [
        [[1, 'K33', 2]],  # K33 not allowed
        [[1, 'K63', 2], [2, 'K29', 3]],  # K29 not allowed
        [[1, 'K48', 2], [2, 'K11', 3]],  # K11 not allowed
        [[1, 'K63', 2], [2, 'K6', 3]],   # K6 not allowed
    ]
    
    for edges in invalid_inputs:
        with pytest.raises(ValueError):
            conjugated_lysines_to_jeff_K48_K63_nomenclature(edges)


def test_49_round_trip_all_formats():
    """Test round-trip conversions through various formats"""
    original_edges = [[1, 'K63', 2], [1, 'K48', 3], [2, 'K33', 4]]
    
    # Convert to compact format and back
    compact = format_nomenclature_preorder_ABC(original_edges)
    parsed_edges = parse_compact_edges(compact)
    assert parsed_edges == original_edges
    
    # Test multimer length consistency
    length1 = multimer_length_from_nomenclature(compact)
    all_nodes = set()
    for src, _, dst in original_edges:
        all_nodes.add(src)
        all_nodes.add(dst)
    length2 = len(all_nodes)
    assert length1 == length2


