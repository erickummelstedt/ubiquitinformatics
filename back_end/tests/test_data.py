five_level_nested_ubiquitin_ = {
    "protein": "1ubq",
    "chain_number": 1,
    "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    "chain_length": 76,
    "branching_sites": [
        {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
        {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
        {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
        {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
        {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
        {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
        {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": {
            "protein": "1ubq",
            "chain_number": 2,
            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "chain_length": 76,
            "branching_sites": [
                {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": {
                    "protein": "1ubq",
                    "chain_number": 3,
                    "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                    "chain_length": 76,
                    "branching_sites": [
                        {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                        {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": {
                            "protein": "1ubq",
                            "chain_number": 4,
                            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                            "chain_length": 76,
                            "branching_sites": [
                                {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                                {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                                {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                                {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                                {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                                {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                                {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                                {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                            ]
                        }},
                        {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": {
                            "protein": "1ubq",
                            "chain_number": 5,
                            "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
                            "chain_length": 76,
                            "branching_sites": [
                                {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
                                {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
                                {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
                                {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                                {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                                {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                                {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                                {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                            ]
                        }},
                        {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
                        {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
                        {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
                        {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
                        {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
                    ]
                }}
            ]
        }},
        {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
    ]
}


# Sample deeply nested ubiquitin dictionary for testing
k48_dimer_ubiquitin = {
    "protein": "1ubq-histag",
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
                                                                                        {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                                                                                        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}]}},
        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}]}

string_k48_dimer_ubiquitin = str(k48_dimer_ubiquitin)

ubiquitin_monomer = {
    "protein": "1ubq",
    "chain_number": 1,
    "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    "chain_length": 76,
    "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                        {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                        {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                        {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                        {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                        {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                        {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}]}

histag_ubiquitin_monomer = {
    "protein": "1ubq",
    "chain_number": 1,
    "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH",
    "chain_length": 83,
    "branching_sites": [{"site_name": "M1","sequence_id": "(M)QIF","children": ""},
                        {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
                        {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
                        {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
                        {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
                        {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
                        {"site_name": "K48","sequence_id": "FAG(K)QLE","children":""}, 
                        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": ""}]}

# Mock helper: minimal working_dictionary
BASE_WORKING_DICT = {
    "protein": "dummy_protein",
    "chain_number": 1,
    "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    "chain_length": 76,
    "branching_sites": [
        {"site_name": "M1", "sequence_id": "(M)QIF", "children": ""},
        {"site_name": "K6", "sequence_id": "IFV(K)TLT", "children": ""},
        {"site_name": "K11", "sequence_id": "LTG(K)TIT", "children": ""},
        {"site_name": "K27", "sequence_id": "ENV(K)AKI", "children": ""},
        {"site_name": "K29", "sequence_id": "VKA(K)IQD", "children": ""},
        {"site_name": "K33", "sequence_id": "IQD(K)EGI", "children": ""},
        {"site_name": "K48", "sequence_id": "FAG(K)QLE", "children": ""},
        {"site_name": "K63", "sequence_id": "NIQ(K)EST", "children": ""}
    ]
}

# Mock helper: minimal context
BASE_CONTEXT = {
    "chain_length_list": [],
    "chain_number_list": [1],
    "free_lysines": [],
    "conjugated_lysines": [],
    "SMAC_lysines": [],
    "ABOC_lysines": [],
    "multimer_string_name": ""
}


ubiquitin_library = {
  "histag_ubi_ubq_1": "{'protein': '1ubq', 'chain_number': 1, 'FASTA_sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH', 'chain_length': 83, 'branching_sites': [{'site_name': 'K48', 'sequence_id': 'FAG(K)QLE', 'children': ''}, {'site_name': 'K11', 'sequence_id': 'LTG(K)TIT', 'children': ''}, {'site_name': 'K63', 'sequence_id': 'NIQ(K)EST', 'children': ''}, {'site_name': 'K6', 'sequence_id': 'IFV(K)TLT', 'children': ''}, {'site_name': 'K27', 'sequence_id': 'ENV(K)AKI', 'children': ''}, {'site_name': 'K29', 'sequence_id': 'VKA(K)IQD', 'children': ''}, {'site_name': 'K33', 'sequence_id': 'IQD(K)EGI', 'children': ''}, {'site_name': 'M1', 'sequence_id': '(M)QIF', 'children': ''}]}",
  "histag_ubi_ubq_1_K48_aboc": "{'protein': '1ubq', 'chain_number': 1, 'FASTA_sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH', 'chain_length': 83, 'branching_sites': [{'site_name': 'K48', 'sequence_id': 'FAG(K)QLE', 'children': 'ABOC'}, {'site_name': 'K11', 'sequence_id': 'LTG(K)TIT', 'children': ''}, {'site_name': 'K63', 'sequence_id': 'NIQ(K)EST', 'children': ''}, {'site_name': 'K6', 'sequence_id': 'IFV(K)TLT', 'children': ''}, {'site_name': 'K27', 'sequence_id': 'ENV(K)AKI', 'children': ''}, {'site_name': 'K29', 'sequence_id': 'VKA(K)IQD', 'children': ''}, {'site_name': 'K33', 'sequence_id': 'IQD(K)EGI', 'children': ''}, {'site_name': 'M1', 'sequence_id': '(M)QIF', 'children': ''}]}",
  "histag_ubi_ubq_1_K63_aboc": "{'protein': '1ubq', 'chain_number': 1, 'FASTA_sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGDHHHHHH', 'chain_length': 83, 'branching_sites': [{'site_name': 'K48', 'sequence_id': 'FAG(K)QLE', 'children': ''}, {'site_name': 'K11', 'sequence_id': 'LTG(K)TIT', 'children': ''}, {'site_name': 'K63', 'sequence_id': 'NIQ(K)EST', 'children': 'ABOC'}, {'site_name': 'K6', 'sequence_id': 'IFV(K)TLT', 'children': ''}, {'site_name': 'K27', 'sequence_id': 'ENV(K)AKI', 'children': ''}, {'site_name': 'K29', 'sequence_id': 'VKA(K)IQD', 'children': ''}, {'site_name': 'K33', 'sequence_id': 'IQD(K)EGI', 'children': ''}, {'site_name': 'M1', 'sequence_id': '(M)QIF', 'children': ''}]}",
  "ubi_ubq_1": "{'protein': '1ubq', 'chain_number': 1, 'FASTA_sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG', 'chain_length': 76, 'branching_sites': [{'site_name': 'K48', 'sequence_id': 'FAG(K)QLE', 'children': ''}, {'site_name': 'K11', 'sequence_id': 'LTG(K)TIT', 'children': ''}, {'site_name': 'K63', 'sequence_id': 'NIQ(K)EST', 'children': ''}, {'site_name': 'K6', 'sequence_id': 'IFV(K)TLT', 'children': ''}, {'site_name': 'K27', 'sequence_id': 'ENV(K)AKI', 'children': ''}, {'site_name': 'K29', 'sequence_id': 'VKA(K)IQD', 'children': ''}, {'site_name': 'K33', 'sequence_id': 'IQD(K)EGI', 'children': ''}, {'site_name': 'M1', 'sequence_id': '(M)QIF', 'children': ''}]}",
  "ubi_ubq_1_K48_SMAC": "{'protein': '1ubq', 'chain_number': 1, 'FASTA_sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG', 'chain_length': 76, 'branching_sites': [{'site_name': 'K48', 'sequence_id': 'FAG(K)QLE', 'children': 'SMAC'}, {'site_name': 'K11', 'sequence_id': 'LTG(K)TIT', 'children': ''}, {'site_name': 'K63', 'sequence_id': 'NIQ(K)EST', 'children': ''}, {'site_name': 'K6', 'sequence_id': 'IFV(K)TLT', 'children': ''}, {'site_name': 'K27', 'sequence_id': 'ENV(K)AKI', 'children': ''}, {'site_name': 'K29', 'sequence_id': 'VKA(K)IQD', 'children': ''}, {'site_name': 'K33', 'sequence_id': 'IQD(K)EGI', 'children': ''}, {'site_name': 'M1', 'sequence_id': '(M)QIF', 'children': ''}]}",
  "ubi_ubq_1_K63_SMAC": "{'protein': '1ubq', 'chain_number': 1, 'FASTA_sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG', 'chain_length': 76, 'branching_sites': [{'site_name': 'K48', 'sequence_id': 'FAG(K)QLE', 'children': ''}, {'site_name': 'K11', 'sequence_id': 'LTG(K)TIT', 'children': ''}, {'site_name': 'K63', 'sequence_id': 'NIQ(K)EST', 'children': 'SMAC'}, {'site_name': 'K6', 'sequence_id': 'IFV(K)TLT', 'children': ''}, {'site_name': 'K27', 'sequence_id': 'ENV(K)AKI', 'children': ''}, {'site_name': 'K29', 'sequence_id': 'VKA(K)IQD', 'children': ''}, {'site_name': 'K33', 'sequence_id': 'IQD(K)EGI', 'children': ''}, {'site_name': 'M1', 'sequence_id': '(M)QIF', 'children': ''}]}",
  "ubi_ubq_1_K48_SMAC_K63_ABOC": "{'protein': '1ubq', 'chain_number': 1, 'FASTA_sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG', 'chain_length': 76, 'branching_sites': [{'site_name': 'K48', 'sequence_id': 'FAG(K)QLE', 'children': 'SMAC'}, {'site_name': 'K11', 'sequence_id': 'LTG(K)TIT', 'children': ''}, {'site_name': 'K63', 'sequence_id': 'NIQ(K)EST', 'children': 'ABOC'}, {'site_name': 'K6', 'sequence_id': 'IFV(K)TLT', 'children': ''}, {'site_name': 'K27', 'sequence_id': 'ENV(K)AKI', 'children': ''}, {'site_name': 'K29', 'sequence_id': 'VKA(K)IQD', 'children': ''}, {'site_name': 'K33', 'sequence_id': 'IQD(K)EGI', 'children': ''}, {'site_name': 'M1', 'sequence_id': '(M)QIF', 'children': ''}]}",
  "ubi_ubq_1_K48_ABOC_K63_SMAC": "{'protein': '1ubq', 'chain_number': 1, 'FASTA_sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG', 'chain_length': 76, 'branching_sites': [{'site_name': 'K48', 'sequence_id': 'FAG(K)QLE', 'children': 'ABOC'}, {'site_name': 'K11', 'sequence_id': 'LTG(K)TIT', 'children': ''}, {'site_name': 'K63', 'sequence_id': 'NIQ(K)EST', 'children': 'SMAC'}, {'site_name': 'K6', 'sequence_id': 'IFV(K)TLT', 'children': ''}, {'site_name': 'K27', 'sequence_id': 'ENV(K)AKI', 'children': ''}, {'site_name': 'K29', 'sequence_id': 'VKA(K)IQD', 'children': ''}, {'site_name': 'K33', 'sequence_id': 'IQD(K)EGI', 'children': ''}, {'site_name': 'M1', 'sequence_id': '(M)QIF', 'children': ''}]}",
  "ubi_ubq_1_K48_ABOC_K63_ABOC": "{'protein': '1ubq', 'chain_number': 1, 'FASTA_sequence': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG', 'chain_length': 76, 'branching_sites': [{'site_name': 'K48', 'sequence_id': 'FAG(K)QLE', 'children': 'ABOC'}, {'site_name': 'K11', 'sequence_id': 'LTG(K)TIT', 'children': ''}, {'site_name': 'K63', 'sequence_id': 'NIQ(K)EST', 'children': 'ABOC'}, {'site_name': 'K6', 'sequence_id': 'IFV(K)TLT', 'children': ''}, {'site_name': 'K27', 'sequence_id': 'ENV(K)AKI', 'children': ''}, {'site_name': 'K29', 'sequence_id': 'VKA(K)IQD', 'children': ''}, {'site_name': 'K33', 'sequence_id': 'IQD(K)EGI', 'children': ''}, {'site_name': 'M1', 'sequence_id': '(M)QIF', 'children': ''}]}"
}

## set the values of the ubiquitin library
ubi_ubq_1_K48_SMAC = ubiquitin_library['ubi_ubq_1_K48_SMAC']
ubi_ubq_1_K63_SMAC = ubiquitin_library['ubi_ubq_1_K63_SMAC']
ubi_ubq_1_K48_SMAC_K63_ABOC = ubiquitin_library['ubi_ubq_1_K48_SMAC_K63_ABOC']
ubi_ubq_1_K48_ABOC_K63_SMAC = ubiquitin_library['ubi_ubq_1_K48_ABOC_K63_SMAC']
ubi_ubq_1_K48_ABOC_K63_ABOC = ubiquitin_library['ubi_ubq_1_K48_ABOC_K63_ABOC']
histag_ubi_ubq_1 = ubiquitin_library['histag_ubi_ubq_1']
histag_ubi_ubq_1_K48_aboc = ubiquitin_library['histag_ubi_ubq_1_K48_aboc']
histag_ubi_ubq_1_K63_aboc = ubiquitin_library['histag_ubi_ubq_1_K63_aboc']

ubi_donor_list = [ubi_ubq_1_K48_SMAC, ubi_ubq_1_K63_SMAC, ubi_ubq_1_K48_SMAC_K63_ABOC, ubi_ubq_1_K48_ABOC_K63_SMAC, ubi_ubq_1_K48_ABOC_K63_ABOC]
ubi_acceptor_list = [histag_ubi_ubq_1, histag_ubi_ubq_1_K48_aboc, histag_ubi_ubq_1_K63_aboc]

