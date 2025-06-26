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
from tests.test_data import *

# =========================================================
# Run tests to ensure the code is working correctly
# This section runs tests to validate the functionality of the code.
# =========================================================

# =========================================================
# Build reactionm database for polyubiquitins
# This section creates a reaction database for polyubiquitin reactions.
# =========================================================

# Function to save reaction database
def save_reaction_database(
        acceptor_list: list, 
        donor_list: list,
        multimer_size: int = 2
    ):
    
    column_names = get_multimer_column_names(multimer_size)

    reaction_histories = create_reaction_histories(acceptor_list, donor_list, multimer_size)

    df = pd.DataFrame(reaction_histories)
    # Expand each list in column 'A' into its own columns
    ubiquitin_history = pd.DataFrame(df['ubiquitin_history'].to_list(), columns=column_names)
    reaction_history = pd.DataFrame(df['reaction_history'].to_list(), columns=column_names)
    donor_history = pd.DataFrame(df['donor_history'].to_list(), columns=column_names)
    context_history = pd.DataFrame(df['context_history'].to_list(), columns=column_names)

    # Save each expanded DataFrame as a CSV file to a relative folder
    output_dir = project_root / 'back_end' / 'data' / 'reaction_database' / f'multimer_size_{multimer_size}'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save the DataFrames to CSV files
    ubiquitin_history.to_csv(output_dir / "ubiquitin_history.csv", index=False)
    reaction_history.to_csv(output_dir / "reaction_history.csv", index=False)
    donor_history.to_csv(output_dir / "donor_history.csv", index=False)
    context_history.to_csv(output_dir / "context_history.csv", index=False)


# Define the acceptor and donor lists
acceptor_list = [
        histag_ubi_ubq_1,
        histag_ubi_ubq_1_K48_aboc,
        histag_ubi_ubq_1_K63_aboc
        ]
donor_list = [
        ubi_ubq_1_K48_SMAC,
        ubi_ubq_1_K63_SMAC,
        ubi_ubq_1_K48_SMAC_K63_ABOC,
        ubi_ubq_1_K48_ABOC_K63_SMAC,
        ubi_ubq_1_K48_ABOC_K63_ABOC
        ]

# Save reaction database for multimer sizes 2 to 6
for i in range(2,6): 
    save_reaction_database(
        acceptor_list=acceptor_list, 
        donor_list=donor_list, 
        multimer_size=i
    )

# =========================================================
# Build base multimers of polyubiquitin
# This section initializes the multimers and builds them up to size 6.
# =========================================================

# Initialize the multimers
multimers = initialize_multimer_dicts(histag_ubi_ubq_1)

# Build multimers of increasing size
for multimer_size in range(2, 6):
    # Expand the multimer list by adding new ubiquitins
    multimers = defining_json_multimers(multimers, ubi_ubq_1)

    # Remove duplicate entries from the multimer dictionary
    delete_duplicate_multimers(multimers)

    # Convert the dictionary to a DataFrame
    multimers_df = pd.DataFrame(multimers)

    # Create output directory for the current multimer size
    output_dir = (
        project_root
        / 'back_end'
        / 'data'
        / 'reaction_database'
        / f'multimer_size_{multimer_size}'
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save the DataFrame as a CSV file
    output_path = output_dir / "multimers.csv"
    multimers_df.to_csv(output_path, index=False)

# =========================================================
# Build reaction database for polyubiquitins with different acceptors and donors
# This section creates a reaction database for polyubiquitin reactions with different acceptors and donors.
# =========================================================     

indexed_values_tetramers, tetarmer_validation_errors = run_script_pulling_indexes_validataing(4, ubiquitin_library)
indexed_values_pentamers, pentamer_validation_errors = run_script_pulling_indexes_validataing(5, ubiquitin_library)

# =========================================================
# Ensure the the building of multimers and reaction database is correct by looking at tetremers
# and pentamers
# Validate indexed values for tetramers and pentamers
# This section checks that the indexed values for tetramers and pentamers match expected values
# and that there are no validation errors.
# =========================================================
# Reference expected values 
# They were selected by manually checking each output
# =========================================================

test_indexed_values_tetramers = [423, 427, 363, 31, 443, 447, 95, 143, 191, 315, 319, 279, 335, 339]
test_indexed_values_pentamers = [
    2035, 2039, 1975, 1655, 47, 2055, 2059, 1719, 111, 1735, 1739, 1039, 127, 871,
    2143, 2147, 2107, 1075, 2163, 2167, 431, 411, 875, 639, 643, 479, 787, 835,
    1491, 1495, 1475, 1307, 1511, 1515, 1311, 1223, 1271, 1599, 1603, 1563, 1619, 1623
]

# Check that actual indexed values match expected values exactly (including order)
if indexed_values_tetramers != test_indexed_values_tetramers:
    raise ValueError(
        f"Indexed values for tetramers do not match expected values.\n"
        f"Expected: {test_indexed_values_tetramers}\n"
        f"Got:      {indexed_values_tetramers}"
    )

if indexed_values_pentamers != test_indexed_values_pentamers: 
    raise ValueError(
        f"Indexed values for pentamers do not match expected values.\n"
        f"Expected: {test_indexed_values_pentamers}\n"
        f"Got:      {indexed_values_pentamers}"
    )

# Check that there are no validation errors
if tetarmer_validation_errors:
    raise ValueError(
        f"Validation errors found for tetramers: {tetarmer_validation_errors}"
    )       
if pentamer_validation_errors:      
    raise ValueError(
        f"Validation errors found for pentamers: {pentamer_validation_errors}"
    )

# ==========================================================
# Check that all the data cleaning, validation and labeling went well
# This section checks that the data cleaning, validation, and labeling processes were successful.
# ==========================================================

tetramer_mismatched_ubiquitin_history, tetramer_mismatched_combined_database = validate_data(4)
pentamer_mismatched_ubiquitin_history, pentamer_mismatched_combined_database = validate_data(5)

# Raise error if mismatched_ubiquitin_history.empty is not True
if not tetramer_mismatched_ubiquitin_history.empty:
    logging.error("Mismatched tetramer ubiquitin history found:")
    logging.error(tetramer_mismatched_ubiquitin_history)
    raise ValueError("Mismatched tetramer ubiquitin history found.")

if not tetramer_mismatched_combined_database.empty:  
    logging.error("Mismatched tetramer combined database found:")
    logging.error(tetramer_mismatched_combined_database)
    raise ValueError("Mismatched tetramer combined database found.")

# Raise error if mismatched_ubiquitin_history.empty is not True
if not pentamer_mismatched_ubiquitin_history.empty:
    logging.error("Mismatched pentamer ubiquitin history found:")
    logging.error(pentamer_mismatched_ubiquitin_history)
    raise ValueError("Mismatched pentamer ubiquitin history found.")

if not pentamer_mismatched_combined_database.empty:  
    logging.error("Mismatched pentamer combined database found:")
    logging.error(pentamer_mismatched_combined_database)
    raise ValueError("Mismatched pentamer combined database found.")

# ==========================================================
# If all checks pass, print success messages
# =========================================================

print("Indexed values for tetramers match expected values.")
print("Indexed values for pentamers match expected values.")
print("No validation errors found for tetramers and pentamers.")

# =========================================================
# If run_file.py completes without errors, it means all tests passed
# and the reaction database and multimers were built correctly.
# =========================================================

print("All tests passed successfully. Reaction database and multimers built correctly.")
# =========================================================