// Example deeply nested ubiquitin JSON for testing
const k48_dimer_ubiquitin = {
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
    {"site_name": "K48",
      "sequence_id": "FAG(K)QLE","children": {
      "protein": "1ubq",
      "chain_number": 2,
      "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
      "chain_length": 76,
      "branching_sites": [
        {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
        {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
        {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
        {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
        {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
        {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
        {"site_name": "K48","sequence_id": "FAG(K)QLE","children": {
      "protein": "1ubq",
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
        {"site_name": "K48","sequence_id": "FAG(K)QLE","children": "ABOC"},
        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}
      ]
    }},
        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": {
      "protein": "1ubq",
      "chain_number": 4,
      "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
      "chain_length": 76,
      "branching_sites": [
        {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
        {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
        {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
        {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
        {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
        {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
        {"site_name": "K48","sequence_id": "FAG(K)QLE","children": "ABOC"},
        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}
      ]
    }}
      ]
    }},
    {"site_name": "K63","sequence_id": "NIQ(K)EST","children": {
      "protein": "1ubq",
      "chain_number": 5,
      "FASTA_sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
      "chain_length": 76,
      "branching_sites": [
        {"site_name": "M1","sequence_id": "(M)QIF","children": ""},
        {"site_name": "K6","sequence_id": "IFV(K)TLT","children": ""},
        {"site_name": "K11","sequence_id": "LTG(K)TIT","children": ""},
        {"site_name": "K27","sequence_id": "ENV(K)AKI","children": ""},
        {"site_name": "K29","sequence_id": "VKA(K)IQD","children": ""},
        {"site_name": "K33","sequence_id": "IQD(K)EGI","children": ""},
        {"site_name": "K48","sequence_id": "FAG(K)QLE","children": "ABOC"},
        {"site_name": "K63","sequence_id": "NIQ(K)EST","children": "SMAC"}
      ]
    }}
  ]
};

export default k48_dimer_ubiquitin;
