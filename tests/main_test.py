import unittest
import json
import logging
from copy import deepcopy



from src.main import (
    iterate_through_ubiquitin,
    inner_wrapper_iterate_through_ubiquitin,
    process_current_protein,
    process_branch,
    handle_protecting_group,
    log_protein_details,
    find_branching_site
)

'''
Explanation of tests for; 
        •	relabelling_ubiquitin_numbers,
        •	inner_wrapper_relabelling_ubiquitin_numbers,
        •	process_current_protein

    1.	Top-Level Tests:
	    •	Validates proper numbering and chain length assignment for the root level.
	2.	Branch Recursion Tests:
	    •	Checks if recursion applies correct numbering and length assignment to nested branches.
	3.	Structure Preservation:
	    •	Ensures site names and sequences remain unchanged after processing.
	4.	Edge Cases:
	    •	Tests the handling of empty children, no branching sites, and invalid input scenarios.
    
	1.	test_relabel_chain_numbers:
	    •	Ensures the chain_number is updated consistently throughout the hierarchy.
	2.	test_no_duplicate_chain_numbers:
	    •	Collects all chain_number values and verifies they are unique, indicating no reuse.
	3.	test_correct_child_structure:
	    •	Validates that child structures (children) contain expected keys and types.
	4.	test_leaf_nodes:
	    •	Ensures all leaf nodes (children == "") are handled properly, with no unexpected branching.

'''

class TestDeeplyNestedProteinStructure(unittest.TestCase):

    def setUp(self):
        """Set up the deeply nested test data."""
        self.deeply_nested_protein = {
            "protein": "1ubq",
            "chain_number": 1,
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGHHHHHH",
            "chain_length": 82,
            "branching_sites": [
                {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                {"site_name": "K48","sequence_id": "FAG(K)QLE","children": {"protein": "1ubq",
                                                                            "chain_number": 2,
                                                                            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                                                            "chain_length": 76,
                                                                            "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                                                                                {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                                                                                {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                                                                                {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                                                                                {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                                                                                {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                                                                                {"site_name": "K48","sequence_id": "FAG(K)QLE","children": {"protein": "1ubq",
                                                                                                                                                            "chain_number": 3,
                                                                                                                                                            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                                                                                                                                            "chain_length": 76,
                                                                                                                                                            "branching_sites": [
                                                                                                                                                                                {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                                                                                                                                                                {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                                                                                                                                                                {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                                                                                                                                                                {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                                                                                                                                                                {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                                                                                                                                                                {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                                                                                                                                                                {"site_name": "K48","sequence_id": "FAG(K)QLE","children": "SMAC"},
                                                                                                                                                                                {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}
                                                                                                                                                                                ]
                                                                                                                                                            }
                                                                                                },
                                                                                                {"site_name": "K63","sequence_id": "NIQ(K)EST","children": {"protein": "1ubq",
                                                                                                                                                            "chain_number": 4,
                                                                                                                                                            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                                                                                                                                            "chain_length": 76,
                                                                                                                                                            "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                                                                                                                                                                {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                                                                                                                                                                {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                                                                                                                                                                {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                                                                                                                                                                {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                                                                                                                                                                {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                                                                                                                                                                {"site_name": "K48","sequence_id": "FAG(K)QLE","children": "SMAC"},
                                                                                                                                                                                {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}
                                                                                                                                                                                ]
                                                                                                                                                            }
                                                                                                }
                                                                                                ]
                                                                                            }
                                                                                            },
                                                                                            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": {"protein": "1ubq",
                                                                                                                                                        "chain_number": 5,
                                                                                                                                                        "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                                                                                                                                        "chain_length": 76,
                                                                                                                                                        "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                                                                                                                                                                            {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                                                                                                                                                                            {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                                                                                                                                                                            {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                                                                                                                                                                            {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                                                                                                                                                                            {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                                                                                                                                                                            {"site_name": "K48","sequence_id": "FAG(K)QLE","children": "SMAC"},
                                                                                                                                                                            {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}
                                                                                                                                                                            ]
                                                                                                                                                        }
                                                                                            }
            ]
        }

    def test_relabelling_ubiquitin_numbers(self):
        """Test top-level relabeling on a deeply nested structure."""
        updated_dict = relabelling_ubiquitin_numbers(self.deeply_nested_protein)
        self.assertIn("chain_number", updated_dict)
        self.assertEqual(updated_dict["chain_number"], 1)
        self.assertEqual(updated_dict["chain_length"], 82)

        # Check first branch's child chain number
        first_branch = updated_dict["branching_sites"][6]["children"]
        self.assertEqual(first_branch["chain_number"], 2)

        # Check deep nested chain number
        deep_nested_chain = (
            first_branch["branching_sites"][7]["children"]
        )
        self.assertEqual(deep_nested_chain["chain_number"], 4)

    def test_inner_wrapper_relabelling_nested(self):
        """Test the recursive relabeling for a deeply nested structure."""
        chain_number_list = [1]
        chain_length_list = []

        updated_dict = inner_wrapper_relabelling_ubiquitin_numbers(
            self.deeply_nested_protein, chain_number_list, chain_length_list
        )

        # Validate top-level changes
        self.assertEqual(updated_dict["chain_number"], 1)
        self.assertEqual(updated_dict["chain_length"], 82)
        self.assertEqual(chain_number_list[1], 2)  # Chain number incremented
        self.assertEqual(len(chain_length_list), 5)

        # Validate a specific branch chain
        first_branch_child = updated_dict["branching_sites"][6]["children"]
        self.assertEqual(first_branch_child["chain_number"], 2)

    def test_branching_sites_structure(self):
        """Ensure branching sites retain their structure after processing."""
        updated_dict = relabelling_ubiquitin_numbers(self.deeply_nested_protein)

        # Traverse to a specific nested branch
        branch_k63 = updated_dict["branching_sites"][6]["children"]["branching_sites"][7]
        self.assertEqual(branch_k63["site_name"], "K63")
        self.assertEqual(branch_k63["children"]["chain_number"], 4)

        # Ensure unchanged branching site names
        self.assertEqual(updated_dict["branching_sites"][0]["site_name"], "M1")
        self.assertEqual(updated_dict["branching_sites"][2]["site_name"], "K11")

    def test_no_duplicate_chain_numbers(self):
        # Call the function
        updated_data = relabelling_ubiquitin_numbers(self.deeply_nested_protein)

        # Collect all chain numbers
        chain_numbers = []

        def extract_chain_numbers(node):
            chain_numbers.append(node['chain_number'])
            for site in node.get('branching_sites', []):
                if site['children'] and isinstance(site['children'], dict):
                    extract_chain_numbers(site['children'])

        extract_chain_numbers(updated_data)

        # Assert all chain numbers are unique
        self.assertEqual(len(chain_numbers), len(set(chain_numbers)))

    def test_empty_children(self):
        """Validate behavior when children are empty strings."""
        updated_dict = relabelling_ubiquitin_numbers(self.deeply_nested_protein)

        # K11 site with no children
        k11_site = updated_dict["branching_sites"][2]
        self.assertEqual(k11_site["site_name"], "K11")
        self.assertEqual(k11_site["children"], "")

    def test_no_branching_sites(self):
        """Ensure function works with no branching sites."""
        test_dict = {
            "protein": "1ubq",
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "branching_sites": []
        }

        updated_dict = relabelling_ubiquitin_numbers(test_dict)
        self.assertEqual(updated_dict["chain_number"], 1)
        self.assertEqual(updated_dict["chain_length"], 76)
        self.assertEqual(updated_dict["branching_sites"], [])

    def test_invalid_input(self):
        """Check behavior on invalid inputs."""
        with self.assertRaises(TypeError):
            relabelling_ubiquitin_numbers(None)
        with self.assertRaises(KeyError):
            relabelling_ubiquitin_numbers({})

if __name__ == "__main__":
    unittest.main()



'''
    1.	test_conjugation_to_empty_lysine:
	    •	Validates that an unoccupied lysine (e.g., K48) can successfully be conjugated with the provided ubiquitin molecule.
	2.	test_no_conjugation_to_occupied_lysine:
	    •	Ensures that attempting to conjugate to an occupied lysine (e.g., K63 with children = "SMAC") raises a TypeError.
	3.	test_invalid_ubi_molecule_sequence:
	    •	Checks that ubiquitin conjugation does not occur if the FASTA_sequence of the ubiquitin molecule is invalid (does not end with 'RLRGG').
	4.	test_update_with_json_input:
	    •	Verifies that the function can handle a JSON string input for the working_dictionary.
	5.	test_logging_message:
	    •	Confirms that the function generates the correct log messages during the conjugation process.
'''

from src.main import K_residue_ubi_addition  # Adjust import as per your module structure

logging.basicConfig(level=logging.INFO)

class TestKResidueUbiAddition(unittest.TestCase):

    def setUp(self):
        # Setup a working dictionary for tests
        self.working_dictionary = {
            "protein": "1ubq",
            "chain_number": 1,
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "branching_sites": [
                {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                { "site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                { "site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                { "site_name": "K33","sequence_id": "IQD(K)EGI","children": "SMAC"},
                { "site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                { "site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}
            ]
        }

        self.working_dictionary_dimer = {
            "protein": "1ubq",
            "chain_number": 1,
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "branching_sites": [
                {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                {"site_name": "K63",
                 "sequence_id": "NIQ(K)EST",
                 "children": {"protein": "1ubq",
                                 "chain_number": 2,
                                 "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                                 "branching_sites": [{ "site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                                                     { "site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                                                     { "site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                                                     { "site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                                                     { "site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                                                     { "site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                                                     { "site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "SMAC"},
                                                     { "site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}]
                             }
                }
            ]
        }

        self.ubi_molecule_to_add = {
            "protein": "1ubq",
            "chain_number": 2,
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "branching_sites": [
                    { "site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                    { "site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                    { "site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                    { "site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                    { "site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                    { "site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                    { "site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "SMAC"},
                    { "site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}]
        }

    def test_conjugation_to_empty_lysine(self):
        # Add ubiquitin to an unoccupied lysine (K48)
        updated_dict = K_residue_ubi_addition(
            working_dictionary=self.working_dictionary,
            specific_ubi_num=1,
            ubiquitination_sequence="FAG(K)QLE",
            ubi_molecule_to_add=self.ubi_molecule_to_add
        )

        # Assert K48 now has the new ubiquitin as its child
        k48_site = updated_dict['branching_sites'][6]
        self.assertEqual(k48_site['children'], self.ubi_molecule_to_add)

    def test_no_conjugation_to_occupied_lysine(self):
        # Attempt to add ubiquitin to an occupied lysine (K63)
        with self.assertRaises(TypeError):
            K_residue_ubi_addition(
                working_dictionary=self.working_dictionary_dimer,
                specific_ubi_num=1,
                ubiquitination_sequence="NIQ(K)EST",
                ubi_molecule_to_add=self.ubi_molecule_to_add
            )

    def test_no_conjugation_to_protected_lysine(self):
        # Attempt to add ubiquitin to an occupied lysine (K63)
        with self.assertRaises(TypeError):
            K_residue_ubi_addition(
                working_dictionary=self.working_dictionary,
                specific_ubi_num=1,
                ubiquitination_sequence="IQD(K)EGI",
                ubi_molecule_to_add=self.ubi_molecule_to_add
            )

    def test_invalid_ubi_molecule_sequence(self):
        # Provide an invalid ubiquitin molecule (missing 'RLRGG' at the end of the sequence)
        invalid_ubi = {
            "protein": "1ubq",
            "chain_number": 2,
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLR",  # Missing 'GG'
            "branching_sites": []
        }

        updated_dict = K_residue_ubi_addition(
            working_dictionary=self.working_dictionary,
            specific_ubi_num=1,
            ubiquitination_sequence="FAG(K)QLE",
            ubi_molecule_to_add=invalid_ubi
        )

        # Assert K48 remains unchanged (conjugation should not happen)
        k48_site = updated_dict['branching_sites'][0]
        self.assertEqual(k48_site['children'], "")

    def test_update_with_json_input(self):
        # Provide a JSON string as input for working_dictionary
        json_input = json.dumps(self.working_dictionary)

        updated_dict = K_residue_ubi_addition(
            working_dictionary=json_input,
            specific_ubi_num=1,
            ubiquitination_sequence="FAG(K)QLE",
            ubi_molecule_to_add=self.ubi_molecule_to_add
        )

        # Assert K48 now has the new ubiquitin as its child
        k48_site = updated_dict['branching_sites'][6]
        self.assertEqual(k48_site['children'], self.ubi_molecule_to_add)

    def test_logging_message(self):
        # Use a test case to check if the correct log messages are generated
        with self.assertLogs(level='INFO') as log:
            K_residue_ubi_addition(
                working_dictionary=self.working_dictionary,
                specific_ubi_num=1,
                ubiquitination_sequence="FAG(K)QLE",
                ubi_molecule_to_add=self.ubi_molecule_to_add
            )
            # Check if specific log messages are present
            self.assertIn("========== ALERT: CONJUGATION ==========", log.output[0])

if __name__ == '__main__':
    unittest.main()


''' 
    1.	SMAC and ABOC Deprotection:
	    •	Verify that protecting groups are removed correctly based on the type of reaction.
	2.	Global Deprotection:
	    •	Test removing all protecting groups in one step.
	3.	K48 and K63 Reactions:
	    •	Validate that ubiquitin is correctly added at K48 or K63.
	4.	Recursive Ubiquitin Chain:
	    •	Test that nested ubiquitin structures are processed recursively.
	5.	Invalid Input Handling:
	    •	Ensure the function raises errors for invalid JSON or missing keys.
'''

from src.main import ubiquitin_simulation

class TestUbiquitinSimulation(unittest.TestCase):

    def setUp(self):
        """Set up test data."""
        # Example parent dictionary (protein/ubiquitin structure)
        self.parent_dict = {
            "protein": "1ubq",
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "branching_sites": [
                {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "SMAC"},
                {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
            ]
        }

        # Example ubiquitin molecule to add
        self.ubi_molecule = {
            "protein": "1ubq",
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "branching_sites": [
                {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
                {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "SMAC"},
                {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
            ]
        }

    def test_smac_deprotection(self):
        """Test SMAC deprotection reaction."""
        expected_result = deepcopy(self.parent_dict)
        expected_result["branching_sites"][5]["children"] = ""
        expected_result["branching_sites"][6]["children"] = ""

        result = ubiquitin_simulation(
            parent_dictionary=self.parent_dict,
            ubi_molecule_to_add=self.ubi_molecule,
            type_of_reaction="SMAC_deprot"
        )
        self.assertEqual(result["branching_sites"], expected_result["branching_sites"])

    def test_aboc_deprotection(self):
        """Test ABOC deprotection reaction."""
        # Add an ABOC protecting group for testing
        parent_dict_with_aboc = deepcopy(self.parent_dict)
        parent_dict_with_aboc["branching_sites"][7]["children"] = "ABOC"

        expected_result = deepcopy(parent_dict_with_aboc)
        expected_result["branching_sites"][7]["children"] = ""

        result = ubiquitin_simulation(
            parent_dictionary=parent_dict_with_aboc,
            ubi_molecule_to_add=self.ubi_molecule,
            type_of_reaction="ABOC_deprot"
        )
        self.assertEqual(result["branching_sites"], expected_result["branching_sites"])

    def test_global_deprotection(self):
        """Test global deprotection reaction."""
        # Add protecting groups for testing
        parent_dict_with_protecting_groups = deepcopy(self.parent_dict)
        parent_dict_with_protecting_groups["branching_sites"][0]["children"] = "SMAC"
        parent_dict_with_protecting_groups["branching_sites"][1]["children"] = "ABOC"

        expected_result = deepcopy(parent_dict_with_protecting_groups)
        for site in expected_result["branching_sites"]:
            site["children"] = ""

        result = ubiquitin_simulation(
            parent_dictionary=parent_dict_with_protecting_groups,
            ubi_molecule_to_add=self.ubi_molecule,
            type_of_reaction="GLOBAL_deprot"
        )
        self.assertEqual(result["branching_sites"], expected_result["branching_sites"])

    def test_k48_conjugation(self):
        """Test K48 reaction (K48 ubiquitin addition)."""
        expected_result = deepcopy(self.parent_dict)
        expected_result["branching_sites"][6]["children"] = deepcopy(self.ubi_molecule)

        input_dictionary = deepcopy(self.parent_dict)
        input_dictionary["branching_sites"][6]["children"] = ''

        result = ubiquitin_simulation(
            parent_dictionary=input_dictionary,
            ubi_molecule_to_add=self.ubi_molecule,
            type_of_reaction="K48"
        )
        self.assertEqual(
            result["branching_sites"][6]["children"]["FASTA_sequence"],
            expected_result["branching_sites"][6]["children"]["FASTA_sequence"]
        )

    def test_k63_conjugation(self):
        """Test K63 reaction (K63 ubiquitin addition)."""
        expected_result = deepcopy(self.parent_dict)
        expected_result["branching_sites"][7]["children"] = deepcopy(self.ubi_molecule)

        result = ubiquitin_simulation(
            parent_dictionary=self.parent_dict,
            ubi_molecule_to_add=self.ubi_molecule,
            type_of_reaction="K63"
        )
        self.assertEqual(
            result["branching_sites"][7]["children"]["FASTA_sequence"],
            expected_result["branching_sites"][7]["children"]["FASTA_sequence"]
        )

    def test_recursive_ubiquitin_chain(self):
        """Test recursive ubiquitin chain addition."""
        # Add an existing ubiquitin to the parent dictionary
        recursive_parent = deepcopy(self.parent_dict)
        recursive_parent["branching_sites"][7]["children"] = deepcopy(self.ubi_molecule)

        expected_result = deepcopy(recursive_parent)
        expected_result["branching_sites"][7]["children"]["branching_sites"] = [
            {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
            {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
            {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
            {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
            {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
            {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": "SMAC"},
            {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": "SMAC"},
            {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": deepcopy(self.ubi_molecule)}
        ]

        result = ubiquitin_simulation(
            parent_dictionary=recursive_parent,
            ubi_molecule_to_add=self.ubi_molecule,
            type_of_reaction="K63"
        )
        self.assertEqual(
            result["branching_sites"][7]["children"]["branching_sites"][7]["children"]["FASTA_sequence"],
            expected_result["branching_sites"][7]["children"]["branching_sites"][7]["children"]["FASTA_sequence"]
        )

    def test_invalid_input(self):
        """Test invalid input raises appropriate errors."""
        with self.assertRaises(json.JSONDecodeError):
            ubiquitin_simulation(
                parent_dictionary="invalid_json",
                ubi_molecule_to_add=self.ubi_molecule,
                type_of_reaction="K48"
            )

        with self.assertRaises(KeyError):
            ubiquitin_simulation(
                parent_dictionary={"invalid_key": "data"},
                ubi_molecule_to_add=self.ubi_molecule,
                type_of_reaction="K48"
            )

if __name__ == "__main__":
    unittest.main()