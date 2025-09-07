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
from src.main import *

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

def combined_history_from_histories(data_dict, ubiquitin_library):

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
    combined_history = pd.concat(
        [donors_mapped_df, reaction_history_filtered, ubiquitin_history_mapped_df],
        axis=0,
        ignore_index=True
    )

    return combined_history

# =========================================
# Data cleaning fucntion
# =========================================

def data_cleaning(combined_history, confirmation_history):
    """
    Perform data cleaning on the combined DataFrame and original data DataFrame.
    Args:
        combined_history (pd.DataFrame): The combined DataFrame with reaction data.
        confirmation_history (pd.DataFrame): The original data DataFrame to be cleaned.
    Returns:   
        pd.DataFrame: The cleaned combined DataFrame and original data DataFrame.
    """
    # Replace Ube13/Mms2_branching with Ubc13/Mms2
    confirmation_history = confirmation_history.map(lambda x: "Ubc13/Mms2" if x == "Ubc13/Mms2 (branching)" else x)
    # Replace Ube13/Mms2 with Ubc13/Mms2
    confirmation_history = confirmation_history.map(lambda x: "Ubc13/Mms2" if x == "Ube13/Mms2" else x)
    # Replace Fake Wash with FAKE_deprot
    confirmation_history = confirmation_history.map(lambda x: "FAKE_deprot" if x == "Fake_Wash" else x)
    # Replace Ubc13/Mms2 (branching) with Ubc13/Mms2 
    combined_history = combined_history.map(lambda x: "Ubc13/Mms2" if x == "Ubc13/Mms2 (branching)" else x)

    return combined_history, confirmation_history

def validate_dataframes_and_extract_indexes_tetramers(combined_history, confirmation_history):
    """ 
    Validate the combined DataFrame and extract indexes for tetramers.
    Args:
        combined_history (pd.DataFrame): The combined DataFrame with reaction data.
        confirmation_history (pd.DataFrame): The original data DataFrame to be cleaned.
    Returns:
        tuple: A tuple containing a list of indexed values and a list of validation errors.
    """
    # Initiate empty list to hold the indexes
    indexed_values = []
    validation_errors = []

    # Iterate through the first 14 rows of confirmation_history
    for i in range(14):
        try:
            # Get the acceptor from the original data
            initial_acceptor = confirmation_history.iloc[(i*3) + 2, 1]
            
            # Get donors from the original data
            dimer_formation_donor = confirmation_history.iloc[(i*3), 2]
            trimer_formation_donor = confirmation_history.iloc[(i*3), 4]
            tetramer_formation_donor = confirmation_history.iloc[(i*3), 6]

            # Get reactions from the original data
            dimer_formation_reaction = confirmation_history.iloc[(i*3)+1, 2]
            dimer_deprotection_reaction = confirmation_history.iloc[(i*3)+1, 3]
            trimer_formation_reaction = confirmation_history.iloc[(i*3)+1, 4]
            trimer_deprotection_reaction = confirmation_history.iloc[(i*3)+1, 5]
            tetramer_formation_reaction = confirmation_history.iloc[(i*3)+1, 6]

            # Filter combined_history
            current_df = combined_history[
                (
                    (combined_history['dimer_formation'] == dimer_formation_reaction) &
                    (combined_history['dimer_deprotection'] == dimer_deprotection_reaction) &
                    (combined_history['trimer_formation'] == trimer_formation_reaction) &
                    (combined_history['trimer_deprotection'] == trimer_deprotection_reaction) &
                    (combined_history['tetramer_formation'] == tetramer_formation_reaction) &
                    (combined_history['table_origin'] == 'Reactions')
                ) |
                (
                    (combined_history['dimer_formation'] == dimer_formation_donor) &
                    (combined_history['trimer_formation'] == trimer_formation_donor) &
                    (combined_history['tetramer_formation'] == tetramer_formation_donor) &
                    (combined_history['table_origin'] == 'Donors')
                ) |
                (
                    (combined_history['initial_acceptor'] == initial_acceptor) &
                    (combined_history['table_origin'] == 'Acceptor')
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

def validate_dataframes_and_extract_indexes_pentamers(combined_history, confirmation_history):
    """ Validate the combined DataFrame and extract indexes for pentamers.
    Args:
        combined_history (pd.DataFrame): The combined DataFrame with reaction data.
        confirmation_history (pd.DataFrame): The original data DataFrame to be cleaned.
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
            initial_acceptor = confirmation_history.iloc[(i*3) + 2, 1]
            
            # Get donors from the original data
            dimer_formation_donor = confirmation_history.iloc[(i*3), 2]
            trimer_formation_donor = confirmation_history.iloc[(i*3), 4]
            tetramer_formation_donor = confirmation_history.iloc[(i*3), 6]
            pentamer_formation_donor = confirmation_history.iloc[(i*3), 8]

            # Get reactions from the original data
            dimer_formation_reaction = confirmation_history.iloc[(i*3)+1, 2]
            dimer_deprotection_reaction = confirmation_history.iloc[(i*3)+1, 3]
            trimer_formation_reaction = confirmation_history.iloc[(i*3)+1, 4]
            trimer_deprotection_reaction = confirmation_history.iloc[(i*3)+1, 5]
            tetramer_formation_reaction = confirmation_history.iloc[(i*3)+1, 6]
            tetramer_deprotecton_reaction = confirmation_history.iloc[(i*3)+1, 7]
            pentamer_formation_reaction = confirmation_history.iloc[(i*3)+1, 8]
            
            # Create a new row for the combined DataFrame   
            current_df = combined_history[
                (
                    (combined_history['dimer_formation'] == dimer_formation_reaction) & \
                    (combined_history['dimer_deprotection'] == dimer_deprotection_reaction) & \
                    (combined_history['trimer_formation'] == trimer_formation_reaction) & \
                    (combined_history['trimer_deprotection'] == trimer_deprotection_reaction) & \
                    (combined_history['tetramer_formation'] == tetramer_formation_reaction) & \
                    (combined_history['tetramer_deprotection'] == tetramer_deprotecton_reaction) & \
                    (combined_history['pentamer_formation'] == pentamer_formation_reaction) & \
                    (combined_history['table_origin'] == 'Reactions')
                )  | \
                (
                    (combined_history['dimer_formation'] == dimer_formation_donor) & \
                    (combined_history['trimer_formation'] == trimer_formation_donor) & \
                    (combined_history['tetramer_formation'] == tetramer_formation_donor) & \
                    (combined_history['pentamer_formation'] == pentamer_formation_donor) & \
                    (combined_history['table_origin'] == 'Donors')
                )  | \
                (
                    (combined_history['initial_acceptor'] == initial_acceptor) & \
                    (combined_history['table_origin'] == 'Acceptor')
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


def assign_ubiquitin_ids(df, index_values, multimer_size):
    """
    Filters and reorders a DataFrame based on a list of index values from a specified column.

    Parameters:
    ----------
    df : pd.DataFrame
        The input DataFrame to filter and reorder.
    index_values : list
        A list of values to retain and reorder by.

    Returns:
    -------
    pd.DataFrame
        A filtered and ordered DataFrame.
    """
    if 'index' not in df.columns:
        raise ValueError(f"Column '{'index'}' not found in DataFrame.")

    filtered = df[df['index'].isin(index_values)]
    ordered = filtered.set_index('index').loc[index_values].reset_index()
    ordered["multimer_id"] = [f"Ub{multimer_size}_{i}" for i in range(1, len(ordered) + 1)]

    return ordered

def merge_multimer_id(df, multimer_numbering_df, indexed_values):
    """
    Merges a DataFrame with the multimer numbering data on 'final_multimer',
    and adds a 'used in synthesis' column based on indexed values and multimer level.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame to merge.
    multimer_numbering_df : pd.DataFrame
        A DataFrame with 'final_multimer' and 'multimer_id' columns.
    indexed_values : list of int
        List of indices used in synthesis (either for tetramers or pentamers).

    Returns
    -------
    pd.DataFrame
        The merged and annotated DataFrame with reordered columns.
    """
    # Merge on 'final_multimer'
    merged = pd.merge(df, multimer_numbering_df, on='index', how='left')

    # Ensure multimer_id is treated as string
    merged['multimer_id'] = merged['multimer_id'].astype(str)

    # Determine expected multimer number from input list length
    if len(indexed_values) == 14:
        expected_multimer_id = '4'
    elif len(indexed_values) == 42:
        expected_multimer_id = '5'
    else:
        raise ValueError("indexed_values must be length 14 (tetramers) or 42 (pentamers)")

    return merged

def mark_used_in_synthesis_by_index(merged: pd.DataFrame, indexed_values: list) -> pd.DataFrame:
    """
    Marks rows in a DataFrame as 'used in synthesis' based on index membership.

    Parameters:
    ----------
    merged : pd.DataFrame
        DataFrame containing an 'index' column.

    indexed_values : list
        List of index values to mark as used in synthesis.

    Returns:
    -------
    pd.DataFrame
        The DataFrame with an added 'used_in_synthesis' column where
        rows with matching indices are marked with 1 and others with 0.
    """
    # Build indicator DataFrame from indexed_values
    synthesis_flags = pd.DataFrame({
        'index': indexed_values,
        'used_in_synthesis': 1
    })

    # Ensure 'index' column is int in both
    merged['index'] = merged['index'].astype(int)
    synthesis_flags['index'] = synthesis_flags['index'].astype(int)

    # Merge on 'index' and fill NaNs with 0
    merged = pd.merge(merged, synthesis_flags, on='index', how='left')
    merged['used_in_synthesis'] = merged['used_in_synthesis'].fillna(0).astype(int)

    # Reorder columns
    cols = merged.columns.tolist()
    reordered_cols = ['index', 'multimer_id', 'used_in_synthesis'] + [col for col in cols if col not in ['index', 'multimer_id', 'used_in_synthesis']]

    return merged[reordered_cols]

# =========================================
# Function to validate confirmation data against the final validation dataset
# =========================================

def tables_are_the_same(confirmation_df, input_df):
    """
    Validate confirmation data against the validation dataset.
    Args:
        confirmation_df (pd.DataFrame): DataFrame containing confirmation data.
        input_df (pd.DataFrame): DataFrame containing validation data.
    Returns:
        pd.DataFrame: DataFrame with validation results.
    """
    # Ensure both DataFrames have the same columns
    confirmation_df = confirmation_df.reindex(columns=input_df.columns, fill_value=np.nan)
    
    # For exact row match (robust):
    matches = confirmation_df.merge(input_df.drop_duplicates(), how='left', indicator=True)
    confirmation_df['is_validated'] = matches['_merge'] == 'both'

    # Step 2: Filter mismatches (optional)
    mismatches = confirmation_df[~confirmation_df['is_validated']]

    return mismatches

# download CSV files
def validate_data(multimer_size):
    """
    Validate the data by loading the combined database and ubiquitin history.
    Args:
        multimer_size (int): Size of the multimer to validate.
    Returns:
        dict: Dictionary containing the combined database and ubiquitin history.
    """

    input_dir = project_root / 'back_end' / 'data' / 'filtered_reaction_database' / f'multimer_size_{multimer_size}'

    # Load the combined database and ubiquitin history
    combined_database = pd.read_csv(input_dir / 'combined_database.csv', index_col=0)
    ubiquitin_history = pd.read_csv(input_dir / 'ubiquitin_history.csv', index_col=0)

    confirmation_dir = project_root / 'back_end' / 'src' / 'confirmation_data' 

    # Load the confirmation data
    confirmation_ubiquitin_history = pd.read_csv(confirmation_dir / f'multimer_size_{multimer_size}_multimer_database.csv')
    confirmation_combined_database = pd.read_csv(confirmation_dir / f'multimer_size_{multimer_size}_reaction_database.csv')

    # Take only the synthesis reactions
    validation_combined_database = combined_database[combined_database['used_in_synthesis']==1]
    validation_ubiquitin_history = ubiquitin_history[ubiquitin_history['used_in_synthesis']==1]

    # Filter columns to match the confirmation database
    checking_columns_ubiquitin_history = confirmation_ubiquitin_history.columns
    validation_ubiquitin_history = validation_ubiquitin_history[checking_columns_ubiquitin_history]

    # Filter columns to match the confirmation database
    checking_columns_combined_database = confirmation_combined_database.columns
    validation_combined_database = validation_combined_database[checking_columns_combined_database]

    # Validate the confirmation data against the validation dataset
    mismatched_ubiquitin_history = tables_are_the_same(
        confirmation_ubiquitin_history, validation_ubiquitin_history
        )

    # Validate the confirmation combined database against the validation dataset
    mismatched_combined_database = tables_are_the_same(
        confirmation_combined_database, validation_combined_database
        )

    return mismatched_ubiquitin_history, mismatched_combined_database


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
    confirmation_history = pd.read_csv(confirmation_data_dir / f"multimer_size_{multimer_size}_reaction_database.csv")

    # =========================================
    # Perform the global deprotection and filtering by SMAC
    # =========================================

    # These two should always go together
    # Perform the global deprotection and filtering
    filtered_data_dict = global_deprotection_filtering_by_smac(data_dict, ubiquitin_library)

    # Create the directory for filtered data
    filtereed_data_dir = project_root / 'back_end' / 'data' / 'filtered_reaction_database' / f'multimer_size_{multimer_size}'
    filtereed_data_dir.mkdir(parents=True, exist_ok=True)

    # Save the filtered histories to CSV files
    ubiquitin_history_filtered = filtered_data_dict['ubiquitin_history']
    reaction_history_filtered = filtered_data_dict['reaction_history']
    donor_history_filtered = filtered_data_dict['donor_history']   
    context_history_filtered = filtered_data_dict['context_history'] 

    # Copy the DataFrames to for iindexes, multimer_ids and whether used in synthesis
    ubiquitin_history_output = ubiquitin_history_filtered.copy()
    reaction_history_output = reaction_history_filtered.copy()
    donor_history_output = donor_history_filtered.copy()
    context_history_output = context_history_filtered.copy()

    # Reset the index and keep the original index as a column
    ubiquitin_history_output.reset_index(inplace=True)
    reaction_history_output.reset_index(inplace=True)
    donor_history_output.reset_index(inplace=True)
    context_history_output.reset_index(inplace=True)

    # Create combined DataFrame from the filtered histories
    # Combine the filtered histories into a single DataFrame
    combined_history = combined_history_from_histories(filtered_data_dict, ubiquitin_library)

    # Perform data cleaning on the combined DataFrame and original data DataFrame
    combined_history, confirmation_history = data_cleaning(combined_history, confirmation_history)

    # Create a copy of the combined DataFrame for indexes, multimer_ids and whether used in synthesis
    combined_history_output = combined_history.copy()
    # Reset the index and keep the original index as a column
    combined_history_output.reset_index(inplace=True)
    # Drop 'level_0' column if it exists
    if 'level_0' in combined_history_output.columns:
        combined_history_output.drop(columns=['level_0'], inplace=True)

    # =========================================
    # Perform validation of synthesis routes chosen
    # =========================================

    # Validate the DataFrames and extract indexes for pentamers
    if multimer_size == 4:
        indexed_values, validation_errors = validate_dataframes_and_extract_indexes_tetramers(combined_history, confirmation_history)
    elif multimer_size == 5:
        # Validate the DataFrames and extract indexes for pentamers
        indexed_values, validation_errors = validate_dataframes_and_extract_indexes_pentamers(combined_history, confirmation_history)

    # =========================================
    # =========================================
    # Save the the dataframes with indexes, multimer_ids and whether used in synthesis
    # =========================================
    # =========================================

    # Assign ubiquitin IDs to the multimers
    ubiquitin_history_synthesised = assign_ubiquitin_ids(ubiquitin_history_output, indexed_values, multimer_size)

    # Step 1: Extract the multimer_numbering_df
    # Very important step as it contains the multimer_id for each final_multimer
    multimer_numbering_synthesised_df = ubiquitin_history_synthesised[['final_multimer', 'multimer_id']].copy()
    # Merge on 'final_multimer'
    merged = pd.merge(ubiquitin_history_output, multimer_numbering_synthesised_df, on='final_multimer', how='left')
    # Create a DataFrame with 'index' and 'multimer_id'
    multimer_numbering_df = merged[['index', 'multimer_id']]

    # Step 2: Merge multimer_id to all target DataFrames
    ubiquitin_history_output = merge_multimer_id(ubiquitin_history_output, multimer_numbering_df, indexed_values)
    reaction_history_output = merge_multimer_id(reaction_history_output, multimer_numbering_df, indexed_values)
    donor_history_output = merge_multimer_id(donor_history_output, multimer_numbering_df, indexed_values)
    context_history_output = merge_multimer_id(context_history_output, multimer_numbering_df, indexed_values)
    combined_history_output = merge_multimer_id(combined_history_output, multimer_numbering_df, indexed_values)

    # Step 3: Mark all the routes that were used in synthesis
    ubiquitin_history_output = mark_used_in_synthesis_by_index(ubiquitin_history_output, indexed_values)
    reaction_history_output = mark_used_in_synthesis_by_index(reaction_history_output, indexed_values)
    donor_history_output = mark_used_in_synthesis_by_index(donor_history_output, indexed_values)
    context_history_output = mark_used_in_synthesis_by_index(context_history_output, indexed_values)
    combined_history_output = mark_used_in_synthesis_by_index(combined_history_output, indexed_values)
    
    # Save the indexed DataFrames to CSV files
    ubiquitin_history_output.to_csv(filtereed_data_dir / "ubiquitin_history.csv", index=True)
    reaction_history_output.to_csv(filtereed_data_dir / "reaction_history.csv", index=True)
    donor_history_output.to_csv(filtereed_data_dir / "donor_history.csv", index=True)
    context_history_output.to_csv(filtereed_data_dir / "context_history.csv", index=True)
    combined_history_output.to_csv(filtereed_data_dir / "combined_database.csv", index=True)

    return indexed_values, validation_errors
















