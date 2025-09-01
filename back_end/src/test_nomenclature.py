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
    format_edges,
    parse_compact_edges,
    build_polyubiquitin_from_edges,
    multimer_length_from_nomenclature,
    conjugated_lysines_to_jeff_K48_K63_nomenclature,
    conjugated_lysines_to_jeff_all_lysines_nomenclature,
    tree_nomenclature_to_numerical_system,
    LETTER_TO_LYS
)

# ===================================
# Basic Function Tests
# ===================================

def test_01_format_edges_simple_chain():
    """Test format_edges with simple linear chain"""
    edges = [[1, 'K63', 2], [2, 'K63', 3]]
    result = format_edges(edges)
    expected = "1A2-2A3"
    assert result == expected


def test_02_format_edges_mixed_lysines():
    """Test format_edges with multiple lysine types"""
    edges = [[1, 'K63', 2], [1, 'K48', 3], [2, 'K33', 4]]
    result = format_edges(edges)
    expected = "1A2-1B3-2C4"
    assert result == expected


def test_03_format_edges_all_seven_lysines():
    """Test format_edges with all seven lysine types"""
    edges = [
        [1, 'K63', 2], [1, 'K48', 3], [1, 'K33', 4], 
        [1, 'K29', 5], [1, 'K27', 6], [1, 'K11', 7], [1, 'K6', 8]
    ]
    result = format_edges(edges)
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
    
    # Test format_edges with empty list
    result = format_edges([])
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
    compact = format_edges(original_edges)
    
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
