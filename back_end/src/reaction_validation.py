import json 
import logging
import copy
import sys
import ast
import numpy as np
from pathlib import Path
import pandas as pd

# Dynamically get the backend path relative to this file
current_file = Path(__file__).resolve()
project_root = current_file.parents[2]  # Go up to project root
sys.path.insert(0, str(project_root))
local_path = project_root / 'back_end'
sys.path.insert(0, str(local_path))

from src.utils.utils import *
from src.utils.logging_utils import *
from main import *

# =========================================
# Functions for DataFrame Manipulation
# =========================================

def global_deprotection_dual(ubiquitin_history_df, context_history_df):
    """
    Applies GLOBAL_deprot to the last column of df1.
    - The first result (multimer) is stored as 'final_multimer' in df1.
    - The second result (context) is stored as the last column in df2.
    """
    last_col_1 = ubiquitin_history_df.columns[-1]
    last_col_2 = context_history_df.columns[-1]

    results = ubiquitin_history_df[last_col_1].apply(lambda x: ubiquitin_simulation(x, '', "GLOBAL_deprot"))

    # Unpack first and second return values into separate series
    ubiquitin_history_df['final_multimer'] = results.apply(lambda r: r[0])
    context_history_df['final_multimer'] = results.apply(lambda r: r[1])

    # Ensure the final_multimer columns are of type string
    ubiquitin_history_df = ubiquitin_history_df.astype(str)
    context_history_df = context_history_df.astype(str)
    
    return ubiquitin_history_df, context_history_df

# =========================================
# Function to filter reaction histories by number of SMAC lysines
# =========================================

def filter_histories_by_number_of_SMAC(
    ubiquitin_history,
    reaction_history,
    donor_history,
    context_history,
    number_of_smac=0,
):
    """
    Filters reaction histories by the maximum number of ABOC_lysines in a specified column
    of the context_history DataFrame.

    Parameters:
    ----------
    ubiquitin_history : pd.DataFrame
        DataFrame containing ubiquitin JSON history.
    reaction_history : pd.DataFrame
        DataFrame containing reaction history.
    donor_history : pd.DataFrame
        DataFrame containing donor species history.
    context_history : pd.DataFrame
        DataFrame of context dictionaries.
    column_offset : int
        Offset from the end of the DataFrame to select the target column (default -2).

    Returns:
    -------
    tuple of pd.DataFrame:
        Filtered versions of (ubiquitin_history, reaction_history,
        donor_history, context_history) based on max ABOC_lysines.
    """

    def extract_aboc_length(cell_string):
        d = ast.literal_eval(cell_string)
        return len(d['ABOC_lysines'])
    
    def extract_smac_length(cell_string):
        d = ast.literal_eval(cell_string)
        return len(d['SMAC_lysines'])
    
    # Apply ABOC length extraction across context_history
    aboc_num_df = context_history.map(extract_aboc_length)

    # Apply SMAC length extraction across context_history
    smac_num_df = context_history.map(extract_smac_length)

    # Select the column of interest using offset
    final_product_column_aboc = aboc_num_df.iloc[:, -2]
    final_product_column_smac = smac_num_df.iloc[:, -2]

    # Identify rows where ABOC count is maximal
    max_value = final_product_column_aboc.max()

    # If number_of_smac is specified, adjust the max_value accordingly
    number_of_aboc = max_value - number_of_smac
     
    filtered_rows = final_product_column_aboc[(final_product_column_aboc == number_of_aboc) & (final_product_column_smac == number_of_smac)]
    selected_indexes = filtered_rows.index

    # Apply index filtering across all input histories
    return (
        ubiquitin_history.loc[selected_indexes],
        reaction_history.loc[selected_indexes],
        donor_history.loc[selected_indexes],
        context_history.loc[selected_indexes]
    )

# =========================================
# Function to move 'table_origin' to the first column in a DataFrame
# =========================================

def move_column_to_front(df, column_name):
    cols = [column_name] + [col for col in df.columns if col != column_name]
    return df[cols]

# =========================================
# Function to validate the most frequent index in a DataFrame
# =========================================

def validate_most_frequent_index(current_df):
    """
    Validates that:
    - There is exactly one most frequent value in 'index' column
    - That value appears exactly 3 times

    Returns:
        most_frequent_value if valid, else None
        error_message if error, else None
    """
    try:
        value_counts = current_df['index'].value_counts()
        most_frequent_count = value_counts.iloc[0]
        top_values = value_counts[value_counts == most_frequent_count]

        if len(top_values) != 1:
            return None, "Not a unique most frequent value."
        if most_frequent_count != 3:
            return None, f"Most frequent count = {most_frequent_count}, expected 3."

        return top_values.index[0], None

    except KeyError:
        return None, "Column 'index' not found in the DataFrame."
    except Exception as e:
        return None, f"Validation failed: {e}"
    
# =========================================
# Function to perform global deprotection and filtering by SMAC
# =========================================

def global_deprotection_filtering_by_smac(data_dict, ubiquitin_library):
    """
    Create a combined DataFrame from the provided data dictionary and ubiquitin library.
    Args:
        data_dict (dict): Dictionary containing DataFrames for ubiquitin history, reaction history, donor history, and context history.
        ubiquitin_library (dict): Dictionary mapping ubiquitin names to their identifiers.
    Returns:
        pd.DataFrame: A combined DataFrame with the specified structure.
    """

    # Extract the DataFrames from the dictionary
    ubiquitin_history = data_dict['ubiquitin_history']
    reaction_history = data_dict['reaction_history']
    donor_history = data_dict['donor_history']
    context_history = data_dict['context_history']

    # Perform the global deprotection operation
    ubiquitin_history, context_history = global_deprotection_dual(ubiquitin_history, context_history)

    ubiquitin_history_filtered, reaction_history_filtered, donor_history_filtered, context_history_filtered = filter_histories_by_number_of_SMAC(
        ubiquitin_history,
        reaction_history,
        donor_history,
        context_history,
        number_of_smac=0
    )

    # Combine the filtered histories into a dictionary (these are filtered by number of SMAC & ABOC)
    filtered_data_dict = {
        "ubiquitin_history": ubiquitin_history_filtered,
        "reaction_history": reaction_history_filtered,
        "donor_history": donor_history_filtered,
        "context_history": context_history_filtered
    }
    return filtered_data_dict

def combined_df_from_histories(data_dict, ubiquitin_library):

    """ 
    Combine the filtered histories into a single DataFrame with specified structure.
    Args:
        data_dict (dict): Dictionary containing DataFrames for ubiquitin history, reaction history, donor history, and context history.
        ubiquitin_library (dict): Dictionary mapping ubiquitin names to their identifiers.                  
    Returns:            
        pd.DataFrame: A combined DataFrame with the specified structure.
    """

    # Extract the DataFrames from the dictionary
    ubiquitin_history_filtered = data_dict['ubiquitin_history']
    reaction_history_filtered = data_dict['reaction_history']
    donor_history_filtered = data_dict['donor_history'] 
    context_history_filtered = data_dict['context_history']

    # Step 1: Reverse the dictionary
    reversed_ubiquitin_library = {v: k for k, v in ubiquitin_library.items()}

    # Example: if your dataframe is called df
    # You can use .applymap() to apply this mapping to every cell
    donors_mapped_df = donor_history_filtered.map(lambda x: reversed_ubiquitin_library.get(x, x))  # fallback to original if not found

    # map ubiquitin history with reversed library
    ubiquitin_history_mapped_df = ubiquitin_history_filtered.map(lambda x: reversed_ubiquitin_library.get(x, x))  # fallback to original if not found

    # Set all values to NaN except for 'initial_acceptor'
    cols_to_null = [col for col in ubiquitin_history_mapped_df.columns if col != 'initial_acceptor']
    ubiquitin_history_mapped_df[cols_to_null] = np.nan

    # Add the 'table_origin' column
    donors_mapped_df['table_origin'] = 'Donors'
    reaction_history_filtered['table_origin'] = 'Reactions'
    ubiquitin_history_mapped_df['table_origin'] = 'Acceptor'

    # Move 'table_origin' to the front for all DataFrames
    donors_mapped_df = move_column_to_front(donors_mapped_df, 'table_origin')
    reaction_history_filtered = move_column_to_front(reaction_history_filtered, 'table_origin')
    ubiquitin_history_mapped_df = move_column_to_front(ubiquitin_history_mapped_df, 'table_origin')

    # Reset index for all DataFrames
    donors_mapped_df.reset_index(drop=False, inplace=True)
    reaction_history_filtered.reset_index(drop=False, inplace=True)
    ubiquitin_history_mapped_df.reset_index(drop=False, inplace=True)

    # Combine the DataFrames into one
    combined_df = pd.concat(
        [donors_mapped_df, reaction_history_filtered, ubiquitin_history_mapped_df],
        axis=0,
        ignore_index=True
    )

    return combined_df

# =========================================
# Data cleaning fucntion
# =========================================

def data_cleaning(combined_df, original_data_df):
    """
    Perform data cleaning on the combined DataFrame and original data DataFrame.
    Args:
        combined_df (pd.DataFrame): The combined DataFrame with reaction data.
        original_data_df (pd.DataFrame): The original data DataFrame to be cleaned.
    Returns:   
        pd.DataFrame: The cleaned combined DataFrame and original data DataFrame.
    """
    # Replace Ube13/Mms2_branching with Ubc13/Mms2
    original_data_df = original_data_df.map(lambda x: "Ubc13/Mms2" if x == "Ubc13/Mms2 (branching)" else x)
    # Replace Ube13/Mms2 with Ubc13/Mms2
    original_data_df = original_data_df.map(lambda x: "Ubc13/Mms2" if x == "Ube13/Mms2" else x)
    # Replace Fake Wash with FAKE_deprot
    original_data_df = original_data_df.map(lambda x: "FAKE_deprot" if x == "Fake_Wash" else x)
    # Replace Ubc13/Mms2 (branching) with Ubc13/Mms2 
    combined_df = combined_df.map(lambda x: "Ubc13/Mms2" if x == "Ubc13/Mms2 (branching)" else x)

    return combined_df, original_data_df

def validate_dataframes_and_extract_indexes_tetramers(combined_df, original_data_df):
    """ 
    Validate the combined DataFrame and extract indexes for tetramers.
    Args:
        combined_df (pd.DataFrame): The combined DataFrame with reaction data.
        original_data_df (pd.DataFrame): The original data DataFrame to be cleaned.
    Returns:
        tuple: A tuple containing a list of indexed values and a list of validation errors.
    """
    # Initiate empty list to hold the indexes
    indexed_values = []
    validation_errors = []

    # Iterate through the first 14 rows of original_data_df
    for i in range(14):
        try:
            # Get the acceptor from the original data
            initial_acceptor = original_data_df.iloc[(i*3) + 2, 1]
            
            # Get donors from the original data
            dimer_formation_donor = original_data_df.iloc[(i*3), 2]
            trimer_formation_donor = original_data_df.iloc[(i*3), 4]
            tetramer_formation_donor = original_data_df.iloc[(i*3), 6]

            # Get reactions from the original data
            dimer_formation_reaction = original_data_df.iloc[(i*3)+1, 2]
            dimer_deprotection_reaction = original_data_df.iloc[(i*3)+1, 3]
            trimer_formation_reaction = original_data_df.iloc[(i*3)+1, 4]
            trimer_deprotection_reaction = original_data_df.iloc[(i*3)+1, 5]
            tetramer_formation_reaction = original_data_df.iloc[(i*3)+1, 6]

            # Filter combined_df
            current_df = combined_df[
                (
                    (combined_df['dimer_formation'] == dimer_formation_reaction) &
                    (combined_df['dimer_deprotection'] == dimer_deprotection_reaction) &
                    (combined_df['trimer_formation'] == trimer_formation_reaction) &
                    (combined_df['trimer_deprotection'] == trimer_deprotection_reaction) &
                    (combined_df['tetramer_formation'] == tetramer_formation_reaction) &
                    (combined_df['table_origin'] == 'Reactions')
                ) |
                (
                    (combined_df['dimer_formation'] == dimer_formation_donor) &
                    (combined_df['trimer_formation'] == trimer_formation_donor) &
                    (combined_df['tetramer_formation'] == tetramer_formation_donor) &
                    (combined_df['table_origin'] == 'Donors')
                ) |
                (
                    (combined_df['initial_acceptor'] == initial_acceptor) &
                    (combined_df['table_origin'] == 'Acceptor')
                )
            ]

            # Run validation
            most_frequent_value, error = validate_most_frequent_index(current_df)

            if error:
                validation_errors.append((i, error))
            else:
                indexed_values.append(int(most_frequent_value))

        except Exception as e:
            validation_errors.append((i, f"Unexpected error: {e}"))

    # Report any validation issues after all iterations
    if validation_errors:
        logging.warning("\nValidation Errors:")
        for i, msg in validation_errors:
            logging.warning(f" - Row {i}: {msg}")

    return indexed_values, validation_errors


# =========================================
# Pentamer formation validation and index extraction
# ========================================= 

def validate_dataframes_and_extract_indexes_pentamers(combined_df, original_data_df):
    """ Validate the combined DataFrame and extract indexes for pentamers.
    Args:
        combined_df (pd.DataFrame): The combined DataFrame with reaction data.
        original_data_df (pd.DataFrame): The original data DataFrame to be cleaned.
    Returns:
        tuple: A tuple containing a list of indexed values and a list of validation errors.
    """


    # Initiate empty list to hold the indexes
    indexed_values = []
    validation_errors = []

    # Take 2 rows from the original data
    for i in range(42): 
        try:
            # Get the acceptor from the original data
            initial_acceptor = original_data_df.iloc[(i*3) + 2, 1]
            
            # Get donors from the original data
            dimer_formation_donor = original_data_df.iloc[(i*3), 2]
            trimer_formation_donor = original_data_df.iloc[(i*3), 4]
            tetramer_formation_donor = original_data_df.iloc[(i*3), 6]
            pentamer_formation_donor = original_data_df.iloc[(i*3), 8]

            # Get reactions from the original data
            dimer_formation_reaction = original_data_df.iloc[(i*3)+1, 2]
            dimer_deprotection_reaction = original_data_df.iloc[(i*3)+1, 3]
            trimer_formation_reaction = original_data_df.iloc[(i*3)+1, 4]
            trimer_deprotection_reaction = original_data_df.iloc[(i*3)+1, 5]
            tetramer_formation_reaction = original_data_df.iloc[(i*3)+1, 6]
            tetramer_deprotecton_reaction = original_data_df.iloc[(i*3)+1, 7]
            pentamer_formation_reaction = original_data_df.iloc[(i*3)+1, 8]
            
            # Create a new row for the combined DataFrame   
            current_df = combined_df[
                (
                    (combined_df['dimer_formation'] == dimer_formation_reaction) & \
                    (combined_df['dimer_deprotection'] == dimer_deprotection_reaction) & \
                    (combined_df['trimer_formation'] == trimer_formation_reaction) & \
                    (combined_df['trimer_deprotection'] == trimer_deprotection_reaction) & \
                    (combined_df['tetramer_formation'] == tetramer_formation_reaction) & \
                    (combined_df['tetramer_deprotection'] == tetramer_deprotecton_reaction) & \
                    (combined_df['pentamer_formation'] == pentamer_formation_reaction) & \
                    (combined_df['table_origin'] == 'Reactions')
                )  | \
                (
                    (combined_df['dimer_formation'] == dimer_formation_donor) & \
                    (combined_df['trimer_formation'] == trimer_formation_donor) & \
                    (combined_df['tetramer_formation'] == tetramer_formation_donor) & \
                    (combined_df['pentamer_formation'] == pentamer_formation_donor) & \
                    (combined_df['table_origin'] == 'Donors')
                )  | \
                (
                    (combined_df['initial_acceptor'] == initial_acceptor) & \
                    (combined_df['table_origin'] == 'Acceptor')
                )
                ]

            # Run validation
            most_frequent_value, error = validate_most_frequent_index(current_df)

            if error:
                validation_errors.append((i, error))
            else:
                indexed_values.append(int(most_frequent_value))

        except Exception as e:
            validation_errors.append((i, f"Unexpected error: {e}"))

    # Report any validation issues after all iterations
    if validation_errors:
        logging.warning("\nValidation Errors:")
        for i, msg in validation_errors:
            logging.warning(f" - Row {i}: {msg}")

    return indexed_values, validation_errors

# Parent function to run the script for pulling indexes and validating the data
def run_script_pulling_indexes_validataing(multimer_size, ubiquitin_library):
    """ 
    This function runs the script to pull indexes and validate the data for a given multimer size.
    It reads the necessary CSV files, performs global deprotection and filtering, cleans the data
    and validates the synthesis routes chosen, and extracts indexes for tetramers or pentamers
    based on the multimer size.
    """

    input_dir = project_root / 'back_end' / 'data' / 'reaction_database' / f'multimer_size_{multimer_size}'

    # Create the dictionary to hold the data
    data_dict = {
        "ubiquitin_history": [],
        "reaction_history": [],
        "donor_history": [],
        "context_history": []
    }
    # Read the CSV files into DataFrames
    data_dict['ubiquitin_history'] = pd.read_csv(input_dir / "ubiquitin_history.csv")
    data_dict['reaction_history'] = pd.read_csv(input_dir / "reaction_history.csv")
    data_dict['donor_history'] = pd.read_csv(input_dir / "donor_history.csv")
    data_dict['context_history'] = pd.read_csv(input_dir / "context_history.csv")

    # open back_end/src/original_data/reaction_summeries/1mer__to_4_reaction_summary.csv
    confirmation_data_dir = project_root / 'back_end' / 'src' / 'confirmation_data' 
    original_data_df = pd.read_csv(confirmation_data_dir / f"multimer_size_{multimer_size}_reaction_database.csv")

    # =========================================
    # Perform the global deprotection and filtering by SMAC
    # =========================================

    # These two should always go together
    # Perform the global deprotection and filtering
    filtered_data_dict = global_deprotection_filtering_by_smac(data_dict, ubiquitin_library)

    # Combine the filtered histories into a single DataFrame
    combined_df = combined_df_from_histories(filtered_data_dict, ubiquitin_library)

    # =========================================
    # Perform data cleaning and validation of synthesis routes chosen
    # =========================================

    # Perform data cleaning on the combined DataFrame and original data DataFrame
    combined_df, original_data_df = data_cleaning(combined_df, original_data_df)

    # Validate the DataFrames and extract indexes for pentamers
    if multimer_size == 4:
        indexed_values, validation_errors = validate_dataframes_and_extract_indexes_tetramers(combined_df, original_data_df)
    elif multimer_size == 5:
        # Validate the DataFrames and extract indexes for pentamers
        indexed_values, validation_errors = validate_dataframes_and_extract_indexes_pentamers(combined_df, original_data_df)

    return indexed_values, validation_errors