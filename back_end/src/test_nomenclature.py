from nomenclature import (
    convert_ubiquitin_nomenclature
    )
# Test cases
test_branching_1 = 'A1B1B2'  # Simple fork: A1 branches to both B1 and B2
test_branching_2 = 'A1B1B2C1C2'  # Double fork: A1→B1,B2 and B1→C1,C2
test_branching_3 = 'A1B1B2C3C4'  # Fork then branch: A1→B1,B2 and B2→C3,C4
test_branching_4 = 'A1B1B2C1C2C3C4'  # Full double fork: both B1 and B2 branch
test_branching_5 = 'A1B1C1C2D1D2'  # B1 branches, then both C positions branch
test_branching_6 = 'A1B1C1D1D2'  # Linear then fork: A1→B1→C1→D1,D2
test_branching_7 = 'A1B2C3C4D7D8'  # Linear then fork: A1→B2→C3,C4→D7,D8
test_branching_8 = 'A1B1B2C2C3D4D5D6'  # Complex: A1→B1,B2→C2,C3→D4,D5,D6
test_branching_9 = 'A1B1C1C2D2D3D4'  # B1→C1,C2 then C2→D2,D3,D4 (invalid)
test_branching_10 = 'A1B1B2C1C3D2D6'  # Mixed: A1→B1,B2→C1,C3→D2,D6
test_branching_11 = 'A1B2C4D7D8'  # Linear then terminal fork
test_branching_12 = 'A1B1B2C1C2C3C4D1D2D3D4D5D6D7D8'  # Maximum branching
test_branching_13 = 'A1B1C3D6'  # Invalid: B1 cannot bind to C3
test_branching_14 = 'A1B1B2C1C4D2D8'  # Mixed valid and invalid branches
test_branching_15 = 'A1B1B2C2C3D5D6D7'  # Complex branching with some invalid
test_branching_16 = 'A1B1B2C3D5'


test_inputs = [test_branching_1, test_branching_2, test_branching_3, test_branching_4, 
               test_branching_5, test_branching_6, test_branching_7, test_branching_8,
               test_branching_9, test_branching_10, test_branching_11, test_branching_12,
               test_branching_13, test_branching_14, test_branching_15, test_branching_16]

print("=== Ubiquitin Nomenclature Converter ===\n")

for i, test_input in enumerate(test_inputs, 1):
    result = convert_ubiquitin_nomenclature(test_input)
    print(f"Input {i}: {result['input']}")
    
    if result['success']:
        print(f"Result: {result['result']}")
    else:
        print(f"Error: {result['result']}")
    print()
