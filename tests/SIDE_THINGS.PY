# --- Tests for Recursive Validation ---

# Use a base valid nested structure for modification
base_nested_ubiquitin = {
    "protein": "1ubq-root",
    "chain_number": 1, # Will be overwritten by iterate_through_ubiquitin
    "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
    "chain_length": 83, # Will be overwritten
    "branching_sites": [
        {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
        {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
        {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
        {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
        {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
        {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
        {"site_name": "K48","sequence_id": "FAG(K)QLE","children": { # NESTED UBIQUITIN STARTS HERE
            "protein": "1ubq-child",
            "chain_number": 99, # Invalid, will be overwritten if validation passes
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "chain_length": 99, # Invalid, will be overwritten
            "branching_sites": [
                {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""},
                {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"} # Protecting group in child
            ]
        }}, # NESTED UBIQUITIN ENDS HERE
        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"} # Protecting group in root
    ]
}

def test_recursive_validation_missing_nested_key():
    """
    Test iterate_through_ubiquitin raises KeyError if a required key is missing
    in a *nested* ubiquitin dictionary.
    """
    invalid_nested_data = deepcopy(base_nested_ubiquitin)
    # Remove a required key from the child dictionary
    del invalid_nested_data["branching_sites"][6]["children"]["FASTA_sequence"]

    with pytest.raises(KeyError, match=r"Missing required keys:.*'FASTA_sequence'.*Allowed keys:"):
        iterate_through_ubiquitin(invalid_nested_data)

def test_recursive_validation_invalid_nested_key():
    """
    Test iterate_through_ubiquitin raises KeyError if an invalid key exists
    in a *nested* ubiquitin dictionary.
    """
    invalid_nested_data = deepcopy(base_nested_ubiquitin)
    # Add an invalid key to the child dictionary
    invalid_nested_data["branching_sites"][6]["children"]["invalid_extra_key"] = "some_value"

    with pytest.raises(KeyError, match=r"Invalid keys found:.*'invalid_extra_key'.*Allowed keys:"):
        iterate_through_ubiquitin(invalid_nested_data)

def test_recursive_validation_invalid_nested_site_name():
    """
    Test iterate_through_ubiquitin raises AssertionError if an invalid site_name
    exists in the branching_sites of a *nested* ubiquitin dictionary.
    """
    invalid_nested_data = deepcopy(base_nested_ubiquitin)
    # Introduce an invalid site name in the child's branching sites
    invalid_nested_data["branching_sites"][6]["children"]["branching_sites"].append(
        {"site_name": "K999", "sequence_id": "XXX(K)XXX", "children": ""}
    )

    # The error message includes the chain number assigned during processing
    with pytest.raises(AssertionError, match=r"Invalid sites found: \['K999'\] in Ubiquitin 2"):
        iterate_through_ubiquitin(invalid_nested_data)

def test_recursive_validation_missing_nested_site_name():
    """
    Test iterate_through_ubiquitin raises AssertionError if a required site_name
    is missing from the branching_sites of a *nested* ubiquitin dictionary.
    """
    invalid_nested_data = deepcopy(base_nested_ubiquitin)
    # Remove a required site (K48) from the child's branching sites
    child_sites = invalid_nested_data["branching_sites"][6]["children"]["branching_sites"]
    invalid_nested_data["branching_sites"][6]["children"]["branching_sites"] = [
        site for site in child_sites if site["site_name"] != "K48"
    ]

    # The error message includes the chain number assigned during processing
    with pytest.raises(AssertionError, match=r"Missing required sites: \['K48'\] in Ubiquitin 2"):
        iterate_through_ubiquitin(invalid_nested_data)

def test_recursive_validation_invalid_nested_sequence_id():
    """
    Test iterate_through_ubiquitin raises AssertionError if an invalid sequence_id
    exists in the branching_sites of a *nested* ubiquitin dictionary.
    """
    invalid_nested_data = deepcopy(base_nested_ubiquitin)
    # Introduce an invalid sequence id in the child's branching sites
    invalid_nested_data["branching_sites"][6]["children"]["branching_sites"].append(
        {"site_name": "K6", "sequence_id": "INVALID(SEQ)", "children": ""} # Duplicate K6 site_name but invalid sequence
    )
    # Need to remove original K6 to avoid duplicate site name error first
    child_sites = invalid_nested_data["branching_sites"][6]["children"]["branching_sites"]
    invalid_nested_data["branching_sites"][6]["children"]["branching_sites"] = [
        site for site in child_sites if site["sequence_id"] != "IFV(K)TLT"
    ]

    # The error message includes the chain number assigned during processing
    with pytest.raises(AssertionError, match=r"Invalid sequences found: \['INVALID\(SEQ\)'\] in Ubiquitin 2"):
         iterate_through_ubiquitin(invalid_nested_data)

def test_recursive_validation_missing_nested_sequence_id():
    """
    Test iterate_through_ubiquitin raises AssertionError if a required sequence_id
    is missing from the branching_sites of a *nested* ubiquitin dictionary.
    """
    invalid_nested_data = deepcopy(base_nested_ubiquitin)
    # Remove a required sequence (FAG(K)QLE) from the child's branching sites
    child_sites = invalid_nested_data["branching_sites"][6]["children"]["branching_sites"]
    invalid_nested_data["branching_sites"][6]["children"]["branching_sites"] = [
        site for site in child_sites if site["sequence_id"] != "FAG(K)QLE"
    ]

    # The error message includes the chain number assigned during processing
    with pytest.raises(AssertionError, match=r"Missing required sequences: \['FAG\(K\)QLE'\] in Ubiquitin 2"):
        iterate_through_ubiquitin(invalid_nested_data)


def test_recursive_validation_nested_site_sequence_mismatch():
    """
    Test iterate_through_ubiquitin raises AssertionError if a site_name and sequence_id
    mismatch occurs in the branching_sites of a *nested* ubiquitin dictionary.
    """
    invalid_nested_data = deepcopy(base_nested_ubiquitin)
    # Swap sequence IDs for K48 and K63 in the child dictionary
    child_sites = invalid_nested_data["branching_sites"][6]["children"]["branching_sites"]
    k48_site = next(site for site in child_sites if site["site_name"] == "K48")
    k63_site = next(site for site in child_sites if site["site_name"] == "K63")
    k48_site["sequence_id"], k63_site["sequence_id"] = k63_site["sequence_id"], k48_site["sequence_id"]

    with pytest.raises(AssertionError) as exc_info:
        iterate_through_ubiquitin(invalid_nested_data)

    # Check for specific mismatch messages within the assertion error
    error_msg = str(exc_info.value)
    assert "site_name: K48, does not correspond with the sequence_id: NIQ(K)EST" in error_msg
    assert "site_name: K63, does not correspond with the sequence_id: FAG(K)QLE" in error_msg
    # The error message includes the chain number assigned during processing
    assert "in Ubiquitin 2" in error_msg # Ensure it's from the nested validation


def test_recursive_validation_invalid_nested_children_format_int():
    """
    Test iterate_through_ubiquitin raises AssertionError if a children field
    has an invalid format (int) in the branching_sites of a *nested* ubiquitin dictionary.
    """
    invalid_nested_data = deepcopy(base_nested_ubiquitin)
    # Set an invalid children format (integer) in the child's branching sites
    invalid_nested_data["branching_sites"][6]["children"]["branching_sites"][0]["children"] = 12345 # M1 site

    # The error message includes the chain number assigned during processing
    with pytest.raises(AssertionError, match=r"Invalid children format: 12345 in Ubiquitin 2"):
        iterate_through_ubiquitin(invalid_nested_data)

def test_recursive_validation_invalid_nested_children_format_string():
    """
    Test iterate_through_ubiquitin raises AssertionError if a children field
    has an invalid string format in the branching_sites of a *nested* ubiquitin dictionary.
    """
    invalid_nested_data = deepcopy(base_nested_ubiquitin)
    # Set an invalid children format (invalid string) in the child's branching sites
    invalid_nested_data["branching_sites"][6]["children"]["branching_sites"][1]["children"] = "INVALID_STRING" # K6 site

    # The error message includes the chain number assigned during processing
    with pytest.raises(AssertionError, match=r"Invalid children format: INVALID_STRING in Ubiquitin 2"):
        iterate_through_ubiquitin(invalid_nested_data)

# --- End of Recursive Validation Tests ---


# --- Additional Edge Cases for iterate_through_ubiquitin ---

def test_iterate_empty_branching_sites_list():
    """
    Test iterate_through_ubiquitin with a valid structure but an empty branching_sites list.
    """
    test_dict = {
        "protein": "1ubq-test",
        "chain_number": 99, # Should be overwritten
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "chain_length": 99, # Should be overwritten
        "branching_sites": [] # Empty list
    }
    updated_dict, context = iterate_through_ubiquitin(test_dict)
    assert updated_dict["chain_number"] == 1
    assert updated_dict["chain_length"] == 76
    assert updated_dict["branching_sites"] == []
    assert context["multimer_string_name"] == "GG-1ubq-test-1-()"
    assert context["free_lysines"] == [] # No K48/K63 in empty list
    assert context["conjugated_lysines"] == []
    assert context["SMAC_lysines"] == []
    assert context["ABOC_lysines"] == []

def test_iterate_malformed_child_missing_branching_sites():
    """
    Test iterate_through_ubiquitin where a nested child dictionary is missing the 'branching_sites' key.
    This should raise a KeyError during the recursive call.
    """
    malformed_child_data = deepcopy(base_nested_ubiquitin)
    # Remove 'branching_sites' from the child dictionary
    del malformed_child_data["branching_sites"][6]["children"]["branching_sites"]

    # Expecting KeyError because validate_protein_keys runs first in the recursive call
    with pytest.raises(KeyError, match=r"Missing required keys:.*'branching_sites'.*Allowed keys:"):
        iterate_through_ubiquitin(malformed_child_data)

def test_iterate_multiple_sibling_branches_with_children():
    """
    Test iterate_through_ubiquitin with a structure where multiple branches at the same level have children.
    Ensures correct numbering and context accumulation.
    """
    multi_branch_data = {
        "protein": "root",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children": deepcopy(ubiquitin_monomer)}, # Child 1 (becomes chain 2)
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": deepcopy(ubiquitin_monomer)}  # Child 2 (becomes chain 3)
        ]
    }
    # Modify children slightly for easier identification
    multi_branch_data["branching_sites"][6]["children"]["protein"] = "child-k48"
    multi_branch_data["branching_sites"][7]["children"]["protein"] = "child-k63"
    # Add a protecting group to one child
    multi_branch_data["branching_sites"][7]["children"]["branching_sites"][7]["children"] = "SMAC" # K63 site in child-k63

    updated_dict, context = iterate_through_ubiquitin(multi_branch_data)

    # Check root
    assert updated_dict["chain_number"] == 1
    assert updated_dict["chain_length"] == 83

    # Check first child (K48 branch)
    child1 = updated_dict["branching_sites"][6]["children"]
    assert child1["protein"] == "child-k48"
    assert child1["chain_number"] == 2
    assert child1["chain_length"] == 76

    # Check second child (K63 branch)
    child2 = updated_dict["branching_sites"][7]["children"]
    assert child2["protein"] == "child-k63"
    assert child2["chain_number"] == 3
    assert child2["chain_length"] == 76

    # Check context
    assert context["chain_number_list"] == [1, 2, 3, 4] # Root, Child1, Child2, Next available
    assert context["chain_length_list"] == [83, 76, 76]
    assert sorted(context["conjugated_lysines"]) == sorted([[1, "K48"], [1, "K63"]])
    # Child 1 has free K48 and K63
    # Child 2 has free K48 and SMAC on K63
    assert sorted(context["free_lysines"]) == sorted([[2, "K48"], [2, "K63"], [3, "K48"]])
    assert sorted(context["SMAC_lysines"]) == sorted([[3, "K63"]])
    assert context["ABOC_lysines"] == []
    assert context["multimer_string_name"] == "his-GG-root-1-(<K48_GG-child-k48-2-()><K63_GG-child-k63-3-(<K63_SMAC>)>)"


def test_iterate_comprehensive_context_check():
    """
    Test iterate_through_ubiquitin with a complex structure and verify all context fields.
    """
    complex_data = {
        "protein": "root-his",
        "chain_number": 1,
        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH", # len 83
        "chain_length": 83,
        "branching_sites": [
            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": "ABOC"}, # ABOC on root
            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
            {"site_name": "K48","sequence_id": "FAG(K)QLE","children": { # Child 1 (chain 2)
                "protein": "child1",
                "chain_number": 2,
                "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG", # len 76
                "chain_length": 76,
                "branching_sites": [
                    {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                    {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                    {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                    {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                    {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                    {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                    {"site_name": "K48","sequence_id": "FAG(K)QLE","children": { # Child 2 (chain 3)
                        "protein": "child2",
                        "chain_number": 3,
                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG", # len 76
                        "chain_length": 76,
                        "branching_sites": [
                            {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                            {"site_name": "K48","sequence_id": "FAG(K)QLE","children": "SMAC"}, # SMAC on child 2
                            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""} # Free K63 on child 2
                        ]
                    }},
                    {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"} # SMAC on child 1
                ]
            }},
            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""} # Free K63 on root
        ]
    }

    updated_dict, context = iterate_through_ubiquitin(complex_data)

    expected_context = {
        "chain_length_list": [83, 76, 76],
        "chain_number_list": [1, 2, 3, 4],
        "free_lysines": sorted([[1, "K63"], [3, "K63"]]), # K63 on root, K63 on child2
        "conjugated_lysines": sorted([[1, "K48"], [2, "K48"]]), # K48 on root, K48 on child1
        "SMAC_lysines": sorted([[2, "K63"], [3, "K48"]]), # K63 on child1, K48 on child2
        "ABOC_lysines": sorted([[1, "K11"]]), # K11 on root
        "multimer_string_name": "his-GG-root-his-1-(<K11_ABOC><K48_GG-child1-2-(<K48_GG-child2-3-(<K48_SMAC>)><K63_SMAC>)>)"
    }

    # Sort lists within context for comparison
    context["free_lysines"].sort()
    context["conjugated_lysines"].sort()
    context["SMAC_lysines"].sort()
    context["ABOC_lysines"].sort()

    assert context == expected_context

def create_deeply_nested_ubiquitin(depth=4):
    """Helper to create a valid deeply nested ubiquitin structure (linear K48 chain)."""
    if depth < 1:
        return None

    # Start with the root (His-tagged)
    root = deepcopy(histag_ubiquitin_monomer)
    root["protein"] = f"root-depth{depth}"
    current_level = root

    # Add nested children
    for i in range(1, depth):
        child = deepcopy(ubiquitin_monomer)
        child["protein"] = f"child-level{i+1}"
        # Find K48 site in the current level and attach the child
        k48_site = next(site for site in current_level["branching_sites"] if site["site_name"] == "K48")
        k48_site["children"] = child
        current_level = child # Move to the newly added child for the next iteration

    # Optional: Add a protecting group at the deepest level for variety
    deepest_k63_site = next(site for site in current_level["branching_sites"] if site["site_name"] == "K63")
    deepest_k63_site["children"] = "SMAC"

    return root

# Create a base structure 4 levels deep (Root + 3 Children)
base_deeply_nested_ubiquitin = create_deeply_nested_ubiquitin(depth=4)

# --- New Tests for Deep Recursive Validation (Level 4) ---

def test_deep_recursive_validation_missing_key_level_4():
    """
    Test iterate_through_ubiquitin raises KeyError if a required key is missing
    in a deeply nested ubiquitin dictionary (Level 4).
    """
    invalid_deep_data = deepcopy(base_deeply_nested_ubiquitin)
    # Navigate to the deepest child (Level 4)
    child1 = invalid_deep_data["branching_sites"][6]["children"] # Level 2
    child2 = child1["branching_sites"][6]["children"]           # Level 3
    child3 = child2["branching_sites"][6]["children"]           # Level 4 (Target)
    # Remove a required key from the deepest child
    del child3["FASTA_sequence"]

    # Error should mention Ubiquitin 4
    with pytest.raises(KeyError, match=r"Missing required keys:.*'FASTA_sequence'.*Allowed keys:.*"):
         iterate_through_ubiquitin(invalid_deep_data)
    # Verify it's specifically KeyError from validate_protein_keys
    try:
        iterate_through_ubiquitin(invalid_deep_data)
    except KeyError as e:
        assert "Missing required keys" in str(e)
    except Exception as e:
        pytest.fail(f"Expected KeyError, but got {type(e).__name__}: {e}")


def test_deep_recursive_validation_invalid_key_level_4():
    """
    Test iterate_through_ubiquitin raises KeyError if an invalid key exists
    in a deeply nested ubiquitin dictionary (Level 4).
    """
    invalid_deep_data = deepcopy(base_deeply_nested_ubiquitin)
    # Navigate to the deepest child (Level 4)
    child1 = invalid_deep_data["branching_sites"][6]["children"] # Level 2
    child2 = child1["branching_sites"][6]["children"]           # Level 3
    child3 = child2["branching_sites"][6]["children"]           # Level 4 (Target)
    # Add an invalid key to the deepest child
    child3["invalid_deep_key"] = "some_value"

    # Error should mention Ubiquitin 4 (though validate_protein_keys doesn't include chain in msg)
    with pytest.raises(KeyError, match=r"Invalid keys found:.*'invalid_deep_key'.*Allowed keys:.*"):
        iterate_through_ubiquitin(invalid_deep_data)
    # Verify it's specifically KeyError from validate_protein_keys
    try:
        iterate_through_ubiquitin(invalid_deep_data)
    except KeyError as e:
        assert "Invalid keys found" in str(e)
    except Exception as e:
        pytest.fail(f"Expected KeyError, but got {type(e).__name__}: {e}")


def test_deep_recursive_validation_invalid_site_name_level_4():
    """
    Test iterate_through_ubiquitin raises AssertionError if an invalid site_name
    exists in the branching_sites of a deeply nested ubiquitin dictionary (Level 4).
    """
    invalid_deep_data = deepcopy(base_deeply_nested_ubiquitin)
    # Navigate to the deepest child (Level 4)
    child1 = invalid_deep_data["branching_sites"][6]["children"] # Level 2
    child2 = child1["branching_sites"][6]["children"]           # Level 3
    child3 = child2["branching_sites"][6]["children"]           # Level 4 (Target)
    # Introduce an invalid site name in the deepest child's branching sites
    child3["branching_sites"].append(
        {"site_name": "K999", "sequence_id": "XXX(K)XXX", "children": ""}
    )

    # Error should mention Ubiquitin 4
    with pytest.raises(AssertionError, match=r"Invalid sites found: \['K999'\] in Ubiquitin 4"):
        iterate_through_ubiquitin(invalid_deep_data)

def test_deep_recursive_validation_missing_site_name_level_4():
    """
    Test iterate_through_ubiquitin raises AssertionError if a required site_name
    is missing from the branching_sites of a deeply nested ubiquitin dictionary (Level 4).
    """
    invalid_deep_data = deepcopy(base_deeply_nested_ubiquitin)
    # Navigate to the deepest child (Level 4)
    child1 = invalid_deep_data["branching_sites"][6]["children"] # Level 2
    child2 = child1["branching_sites"][6]["children"]           # Level 3
    child3 = child2["branching_sites"][6]["children"]           # Level 4 (Target)
    # Remove a required site (e.g., K11) from the deepest child's branching sites
    child3["branching_sites"] = [
        site for site in child3["branching_sites"] if site["site_name"] != "K11"
    ]

    # Error should mention Ubiquitin 4
    with pytest.raises(AssertionError, match=r"Missing required sites: \['K11'\] in Ubiquitin 4"):
        iterate_through_ubiquitin(invalid_deep_data)

def test_deep_recursive_validation_invalid_sequence_id_level_4():
    """
    Test iterate_through_ubiquitin raises AssertionError if an invalid sequence_id
    exists in the branching_sites of a deeply nested ubiquitin dictionary (Level 4).
    """
    invalid_deep_data = deepcopy(base_deeply_nested_ubiquitin)
    # Navigate to the deepest child (Level 4)
    child1 = invalid_deep_data["branching_sites"][6]["children"] # Level 2
    child2 = child1["branching_sites"][6]["children"]           # Level 3
    child3 = child2["branching_sites"][6]["children"]           # Level 4 (Target)
    # Introduce an invalid sequence id in the deepest child's branching sites
    # Replace K6's sequence
    k6_site = next(site for site in child3["branching_sites"] if site["site_name"] == "K6")
    k6_site["sequence_id"] = "INVALID(SEQ)"

    # Error should mention Ubiquitin 4
    with pytest.raises(AssertionError, match=r"Invalid sequences found: \['INVALID\(SEQ\)'\] in Ubiquitin 4"):
         iterate_through_ubiquitin(invalid_deep_data)

def test_deep_recursive_validation_missing_sequence_id_level_4():
    """
    Test iterate_through_ubiquitin raises AssertionError if a required sequence_id
    is missing from the branching_sites of a deeply nested ubiquitin dictionary (Level 4).
    """
    invalid_deep_data = deepcopy(base_deeply_nested_ubiquitin)
    # Navigate to the deepest child (Level 4)
    child1 = invalid_deep_data["branching_sites"][6]["children"] # Level 2
    child2 = child1["branching_sites"][6]["children"]           # Level 3
    child3 = child2["branching_sites"][6]["children"]           # Level 4 (Target)
    # Remove a required sequence (e.g., K27's sequence) from the deepest child's branching sites
    child3["branching_sites"] = [
        site for site in child3["branching_sites"] if site["sequence_id"] != "ENV(K)AKI"
    ]

    # Error should mention Ubiquitin 4
    with pytest.raises(AssertionError, match=r"Missing required sequences: {'ENV(K)AKI'} in Ubiquitin 4"):
        iterate_through_ubiquitin(invalid_deep_data)


def test_deep_recursive_validation_site_sequence_mismatch_level_4():
    """
    Test iterate_through_ubiquitin raises AssertionError if a site_name and sequence_id
    mismatch occurs in the branching_sites of a deeply nested ubiquitin dictionary (Level 4).
    """
    invalid_deep_data = deepcopy(base_deeply_nested_ubiquitin)
    # Navigate to the deepest child (Level 4)
    child1 = invalid_deep_data["branching_sites"][6]["children"] # Level 2
    child2 = child1["branching_sites"][6]["children"]           # Level 3
    child3 = child2["branching_sites"][6]["children"]           # Level 4 (Target)
    # Swap sequence IDs for K27 and K29 in the deepest child dictionary
    k27_site = next(site for site in child3["branching_sites"] if site["site_name"] == "K27")
    k29_site = next(site for site in child3["branching_sites"] if site["site_name"] == "K29")
    k27_site["sequence_id"], k29_site["sequence_id"] = k29_site["sequence_id"], k27_site["sequence_id"]

    with pytest.raises(AssertionError) as exc_info:
        iterate_through_ubiquitin(invalid_deep_data)

    # Check for specific mismatch messages within the assertion error
    error_msg = str(exc_info.value)
    assert "site_name: K27, does not correspond with the sequence_id: VKA(K)IQD" in error_msg
    assert "site_name: K29, does not correspond with the sequence_id: ENV(K)AKI" in error_msg
    # Error should mention Ubiquitin 4 (validate_branching_sites adds this)
    assert "in Ubiquitin 4" in error_msg # Ensure it's from the nested validation


def test_deep_recursive_validation_invalid_children_format_int_level_4():
    """
    Test iterate_through_ubiquitin raises AssertionError if a children field
    has an invalid format (int) in the branching_sites of a deeply nested ubiquitin dictionary (Level 4).
    """
    invalid_deep_data = deepcopy(base_deeply_nested_ubiquitin)
    # Navigate to the deepest child (Level 4)
    child1 = invalid_deep_data["branching_sites"][6]["children"] # Level 2
    child2 = child1["branching_sites"][6]["children"]           # Level 3
    child3 = child2["branching_sites"][6]["children"]           # Level 4 (Target)
    # Set an invalid children format (integer) in the deepest child's branching sites (e.g., M1 site)
    m1_site = next(site for site in child3["branching_sites"] if site["site_name"] == "M1")
    m1_site["children"] = 98765

    # Error should mention Ubiquitin 4
    with pytest.raises(AssertionError, match=r"Invalid children format: 98765"):
        iterate_through_ubiquitin(invalid_deep_data)

def test_deep_recursive_validation_invalid_children_format_string_level_4():
    """
    Test iterate_through_ubiquitin raises AssertionError if a children field
    has an invalid string format in the branching_sites of a deeply nested ubiquitin dictionary (Level 4).
    """
    invalid_deep_data = deepcopy(base_deeply_nested_ubiquitin)
    # Navigate to the deepest child (Level 4)
    child1 = invalid_deep_data["branching_sites"][6]["children"] # Level 2
    child2 = child1["branching_sites"][6]["children"]           # Level 3
    child3 = child2["branching_sites"][6]["children"]           # Level 4 (Target)
    # Set an invalid children format (invalid string) in the deepest child's branching sites (e.g., K6 site)
    k6_site = next(site for site in child3["branching_sites"] if site["site_name"] == "K6")
    k6_site["children"] = "NOT_VALID_PG"

    # Error should mention Ubiquitin 4
    with pytest.raises(AssertionError, match=r"Invalid children format: NOT_VALID_PG"):
        iterate_through_ubiquitin(invalid_deep_data)

# --- End of Deep Recursive Validation Tests ---