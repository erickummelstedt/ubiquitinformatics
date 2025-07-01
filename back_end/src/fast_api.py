import json 
import logging
import copy
import sys
import ast
import numpy as np
from pathlib import Path
import pandas as pd
import ast 

import matplotlib.pyplot as plt
from matplotlib import font_manager

# Dynamically get the backend path relative to this file
current_file = Path(__file__).resolve()
project_root = current_file.parents[2]  # Go up to project root
sys.path.insert(0, str(project_root))
local_path = project_root / 'back_end'
sys.path.insert(0, str(local_path))


def plot_96wells(figure=1, figure_name = 'Test',colorbar_type= 'PuRd', cdata=None, sdata=None, bdata=None, bcolors=None, bmeans=None, **kwargs):
    # from https://github.com/jaumebonet/RosettaSilentToolbox/blob/master/rstoolbox/plot/experimental.py
    """Plot data of a 96 well plate into an equivalent-shaped plot.

    Allows up to three data sources at the same time comming from experiments performed in a
    `96 wells microplate <https://www.wikiwand.com/en/Microtiter_plate>`_. Two of this data
    sources can have to be numerical, and are represented as color and size, while an extra
    boolean dataset can be represented by the color of the border of each well.

    Plot is based on :func:`~matplotlib.pyplot.scatter`; some graphical properties to control
    the visuals (such as ``cmap``), can be provided through this function.

    :param cdata: Data contentainer to be represented by color coding. Has to contain 8 rows
        and 12 columns, like the 96 well plate. Contains **continuos numerical data**.
    :type cdata: :class:`~pandas.DataFrame`
    :param sdata: Data contentainer to be represented by size coding. Has to contain 8 rows
        and 12 columns, like the 96 well plate. Contains **continuos numerical data**.
    :type sdata: :class:`~pandas.DataFrame`
    :param bdata: Data contentainer to be represented by the edge color. Has to contain 8 rows
        and 12 columns, like the 96 well plate. Contains **boolean data**.
    :type bdata: :class:`~pandas.DataFrame`
    :param bcolors: List of the two colors to identify the differences in the border for binary
        data. It has to be a list of 2 colors only. First color represents :data:`True` instances
        while the second color is :data:`False`. Default is ``black`` and ``green``.
    :type bcolors: :func:`list` of :class:`str`
    :param bmeans: List with the meanings of the boolean condition (for the legend). First color
        represents :data:`True` instances while the second color is :data:`False`. Default is
        ``True`` and ``False``.
    :type bmeans: :func:`list` of :class:`str`

    :return: Union[:class:`~matplotlib.figure.Figure`, :class:`~matplotlib.axes.Axes`]

    :raises:
        :ValueError: If input data is not a :class:`~pandas.DataFrame`.
        :ValueError: If :class:`~pandas.DataFrame` do not has the proper shape.
        :ValueError: If ``bcolors`` of ``bmeans`` are provided with sizes different than 2.

    .. rubric:: Example

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.plot import plot_96wells
           ...: import numpy as np
           ...: import pandas as pd
           ...: import matplotlib.pyplot as plt
           ...: np.random.seed(0)
           ...: df = pd.DataFrame(np.random.randn(8, 12))
           ...: fig, ax = plot_96wells(cdata = df, sdata = -df, bdata = df<0)
           ...: plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)

        @savefig plot_96wells_docs.png width=5in
        In [2]: plt.show()

        In [3]: plt.close()
    """

    # Changes in one of this parameters should change others to ensure size fit.
    fig = plt.figure(figure, figsize=(15, 7))
    fig.suptitle(figure_name, fontsize=16)
    ax  = plt.subplot2grid((1, 1), (0, 0), fig=fig)
    top = 2000
    bot = 50
    sizeOfFont = 15
    ticks_font = font_manager.FontProperties(style='normal', size=sizeOfFont, weight='normal')

    # Fixed: THIS CANNOT CHANGE!
    kwargs['x'] = list(range(1, 13)) * 8
    kwargs['y'] = sorted(list(range(1, 9)) * 12)

    # Modified by the input data.
    kwargs.setdefault('s', [top, ] * len(kwargs['y']))
    kwargs.setdefault('c', 'white')
    kwargs.setdefault('edgecolor', ['black', ] * len(kwargs['y']))
    kwargs.setdefault('linewidths', 1.5)
    kwargs.setdefault('cmap', colorbar_type)

    # Color Data
    if cdata is not None:
        if not isinstance(cdata, pd.DataFrame):
            raise ValueError('Wrong data type')
        if cdata.shape != (8, 12):
            raise ValueError('Wrong data shape')
        kwargs['c'] = cdata.values.flatten()

    # Size Data
    if sdata is not None:
        if not isinstance(sdata, pd.DataFrame):
            raise ValueError('Wrong data type')
        if sdata.shape != (8, 12):
            raise ValueError('Wrong data shape')
        p = sdata.values.flatten()
        p = ((p - np.min(p)) / np.ptp(p)) * (top - bot) + bot
        kwargs['s'] = p

    # Border Data
    if bdata is not None:
        if not isinstance(bdata, pd.DataFrame):
            raise ValueError('Wrong data type')
        if bdata.shape != (8, 12):
            raise ValueError('Wrong data shape')
        if not (bdata.dtypes == bool).all():
            raise ValueError('bdata has to be booleans')
        if bcolors is None:
            bcolors = ['black', 'green']
        if bmeans is None:
            bmeans = ['True', 'False']
        if len(bcolors) < 2:
            raise ValueError('At least to border colors need to be provided')
        if len(bmeans) < 2:
            raise ValueError('At least to binary names need to be provided')
        b = bdata.values.flatten()
        b = [bcolors[0] if _ else bcolors[1] for _ in b]
        kwargs['edgecolor'] = b

    # PLOT
    mesh = ax.scatter(**kwargs)

    # Image aspects
    ax.grid(False)
    ax.set_xticks(range(1, 13))
    ax.xaxis.tick_top()
    ax.set_yticks(range(1, 9))
#    ax.set_yticklabels(string.ascii_uppercase[0:9])
    ax.set_ylim((8.5, 0.48))
    ax.set_aspect(1)
    ax.tick_params(axis='both', which='both', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
    
    ## adding text to 96 well plate
    for text, x, y in zip(cdata.values.flatten(), list(range(1, 13)) * 8, sorted(list(range(1, 9)) * 12)):
        if text == 0: 
            new_text = ''
        else:
            new_text = text 
        if text < 10:
            plt.text(x-0.175, y+0.15, str(new_text), fontsize = 23)
        else:
            plt.text(x-0.275, y+0.15, str(new_text), fontsize = 23)

    return fig, ax


# ============================================
# Functions 
# ============================================

# This function creates an encoded dictionary that maps combinations of donors and enzymes to unique integers.
def create_e_d_encoded_dictionary(list_of_enzymes, list_of_donors):
    """
    Creates a dictionary that encodes the combinations of donors and enzymes.
    
    Parameters:
        list_of_enzymes (list): List of enzyme names.
        list_of_donors (list): List of donor names.
        
    Returns:
        dict: Dictionary mapping (donor, enzyme) tuples to unique integers.
    """
    e_d_encoded_dictionary = {}
    count = 1 
    for i in list_of_enzymes: 
        for j in list_of_donors:
            ## these are the self reactions should not be included
            if ((j == 'ubi_ubq_1_K48_SMAC') & (i=='Ubc13/Mms2')) | ((j == 'ubi_ubq_1_K63_SMAC') & (i=='Ube2K'))|((j == 'ubi_ubq_1_K63_SMAC') & (i=='gp78/Ube2g2')):
                count = count
            else:
                e_d_encoded_dictionary[j, i] = count 
                count = count+1 
    return e_d_encoded_dictionary  

def select_columns_by_keyword(combined_database, keyword):
    """
    Select columns from a DataFrame that contain a specific keyword (case-insensitive).

    Parameters:
        combined_database (pd.DataFrame): The input DataFrame.
        keyword (str): The keyword to search for in column names.

    Returns:
        list: A list of matching column names.
    """
    return [col for col in combined_database.columns if keyword.lower() in col.lower()]


def get_combined_data(combined_database, index_list, table_origin, keyword='formation'):
    """
    Returns formation-related columns for specified indices and table_origin,
    preserving the order of index_list.
    
    Parameters:
        combined_database (pd.DataFrame): The full database to filter from.
        index_list (list of int): List of index values to retrieve.
        table_origin (str): 'Reactions', 'Donors', etc.
        
    Returns:
        pd.DataFrame: Filtered DataFrame with only formation-related columns.
    """
    # Filter by condition
    filtered_rows = combined_database[
        (combined_database['index'].isin(index_list)) &
        (combined_database['table_origin'] == table_origin)
    ]

    # Preserve order of index_list
    filtered_rows_sorted = filtered_rows.set_index('index').loc[index_list].reset_index()

    # Select columns that contain 'formation' in the name
    formation_columns = select_columns_by_keyword(combined_database, keyword)
    
    # Extract only those columns
    return filtered_rows_sorted[formation_columns]

# Modified folding function: always creates two columns, padding with None if necessary
def fold_series(series, chunk_length=8, min_chunks=2):
    total_length = len(series)
    num_chunks = max(int(np.ceil(total_length / chunk_length)), min_chunks)
    padded = series.tolist() + [None] * (num_chunks * chunk_length - total_length)
    folded = [padded[i * chunk_length:(i + 1) * chunk_length] for i in range(num_chunks)]
    return pd.DataFrame(folded).transpose()

# Fold DataFrame function: applies folding to each column and renames columns for clarity
def fold_dataframe(df):
    """
    Apply folding to each column in a DataFrame using a provided fold_series function,
    then rename columns for clarity.

    Parameters:
        df (pd.DataFrame): The input DataFrame.
        fold_series (function): Function to apply to each column (e.g., for encoding).

    Returns:
        pd.DataFrame: The folded and renamed DataFrame.
    """
    # Check for maximum length constraint
    MAX_LENGTH = 16
    if len(df) > MAX_LENGTH:
        raise ValueError(f"DataFrame length {len(df)} exceeds maximum allowed length of {MAX_LENGTH}.")

    # Apply folding to each column
    folded_columns = [fold_series(df[col]) for col in df.columns]

    # Concatenate folded columns
    enzymes_donors_96 = pd.concat(folded_columns, axis=1, keys=df.columns)

    # Rename columns for clarity
    enzymes_donors_96.columns = [f"{col}_{i+1}" for col, i in enzymes_donors_96.columns]

    return enzymes_donors_96


def pad_and_validate_dataframe(input_df, required_columns=12):
    """
    Pad the DataFrame with empty columns if it has fewer than required_columns.
    Raise an error if it exceeds the limit. Convert NaNs to 0 and cast to integers.

    Parameters:
        input_df (pd.DataFrame): The input DataFrame to pad and validate.
        required_columns (int): The required number of columns (default: 12).

    Returns:
        pd.DataFrame: A cleaned and padded DataFrame with exactly required_columns columns.
    """
    current_columns = input_df.shape[1]

    # If fewer than required, add empty (None) columns
    if current_columns < required_columns:
        for i in range(current_columns + 1, required_columns + 1):
            input_df[f"padding_{i}"] = None

    # If more than required, raise an error
    elif current_columns > required_columns:
        raise ValueError(f"Final DataFrame has {current_columns} columns, which exceeds the {required_columns}-column limit.")

    # Convert the DataFrame to integers, filling NaNs with 0
    input_df = input_df.fillna(0).infer_objects().astype(int)

    return input_df

def inner_create_plate_dfs(data_dict, indexed_values, multimer_size=4):
    """    
    Creates a DataFrame with formation data for dimer, trimer, and tetramer levels,
    including encoded values based on donor and reaction combinations.  
    Parameters:
        combined_database (pd.DataFrame): The full database to filter from.
        indexed_values (list of int): List of index values for tetramers.
        e_d_encoded_dictionary (dict): Dictionary mapping (donor, reaction) tuples to
        unique integers.

    Returns:
        dict: A dictionary containing the encoded DataFrame and the encoded dictionary.
        The dictionary will have the following keys:
            - 'enzymes_donors_96': DataFrame with encoded donor and reaction formation
              data for dimer, trimer, and tetramer levels.
            - 'deprots_96': DataFrame with encoded deprotection data for dimer and trimer levels.
            - 'dimer_acceptors_96': DataFrame with encoded dimer acceptor data
              for dimer levels.
            - 'e_d_encoded_dictionary': Dictionary mapping (donor, reaction) tuples to
              unique integers for donor and reaction combinations.
            - 'deprot_encoded_dictionary': Dictionary mapping deprotection types to unique integers.
            - 'dimer_encoded_dictionary': Dictionary mapping dimer acceptor names to unique integers.
    """    
    # ================================
    # Extract the necessary data from the data_dict
    # ================================
    combined_database = data_dict['combined_database']
    context_history = data_dict['context_history']

    # ================================
    # Determine the levels of formation based on multimer size
    # ================================
    if multimer_size == 4:
        levels = ['dimer', 'trimer', 'tetramer']
    elif multimer_size == 5:
        levels = ['dimer', 'trimer', 'tetramer', 'pentamer']

    # ================================
    # Reactants and Donors
    # ================================
    # Inputs:
    # combined_database: The full database to filter from.
    # indexed_values: List of index values for tetramers.
    # e_d_encoded_dictionary: Dictionary mapping (donor, reaction) tuples to unique integers.
    # ================================

    # List of donors and enzymes to encode
    list_of_donors = [
        'ubi_ubq_1_K48_SMAC', 
        'ubi_ubq_1_K63_SMAC', 
        'ubi_ubq_1_K48_SMAC_K63_ABOC', 
        'ubi_ubq_1_K48_ABOC_K63_SMAC', 
        'ubi_ubq_1_K48_ABOC_K63_ABOC'
        ]

    list_of_enzymes = [
        'gp78/Ube2g2', 
        'Ube2K', 
        'Ubc13/Mms2'
        ]
    
    # Create the encoded dictionary
    e_d_encoded_dictionary = create_e_d_encoded_dictionary(list_of_enzymes, list_of_donors)

    # ================================
     # Create an empty DataFrame to hold the results
    enzymes_donors_96 = pd.DataFrame()

    # Filter the combined database for reactions and donors based on indexed values
    filtered_reactions = get_combined_data(combined_database, indexed_values, 'Reactions', 'formation')
    filtered_donors = get_combined_data(combined_database, indexed_values, 'Donors', 'formation')

    # Iterate through each level of formation
    for level in levels:
        if level == 'dimer':
            continue  # Skip this iteration
        else:
            # Create columns for donor and reaction formation data
            enzymes_donors_96[f'{level}_donor'] = filtered_donors[f'{level}_formation']
            enzymes_donors_96[f'{level}_reaction'] = filtered_reactions[f'{level}_formation']

            # Encode the donor and reaction combinations using the encoded dictionary
            # e_d_encoded mean enzyme donor encoded
            enzymes_donors_96[f'{level}_e_d_encoded'] = enzymes_donors_96.apply(
                lambda row: e_d_encoded_dictionary.get((row[f'{level}_donor'], row[f'{level}_reaction']), None), axis=1
            )

    # Pull only the encoded columns
    e_d_encoded_columns = select_columns_by_keyword(enzymes_donors_96, '_e_d_encoded')

    # Filter the DataFrame to keep only the encoded columns
    enzymes_donors_96 = enzymes_donors_96[e_d_encoded_columns]

    # Fold the DataFrame to ensure each column has a maximum of 8 values
    # Final DataFrame should have 12 columns, padded if necessary with encoded values
    enzymes_donors_96 = fold_dataframe(enzymes_donors_96)
    
    # ================================
    # Functions for Deprotections
    # ================================
    # Inputs:
    # filtered_reactions: The filtered reactions DataFrame containing deprotection information.
    # indexed_values: List of index values for tetramers.
    # levels: List of levels for deprotection (e.g., dimer, trimer, tetramer).
    # ================================
    
    # Create the encoded dictionary
    deprot_encoded_dictionary = {
        'SMAC_deprot': 1,
        'FAKE_deprot': 2
        }
    
    # Create an empty DataFrame to hold the results
    deprots_96 = pd.DataFrame()

    filtered_deprots = get_combined_data(combined_database, indexed_values, 'Reactions', 'deprotection')

    for level in levels[:-1]:  # Exclude 'final' level as the reactions end with a formation
        # Create columns for donor and reaction formation data
        deprots_96[f'{level}_deprotection'] = filtered_deprots[f'{level}_deprotection']

        # Encode the donor and reaction combinations using the encoded dictionary
        # e_d_encoded mean enzyme donor encoded
        deprots_96[f'{level}_deprot_encoded'] = deprots_96.apply(
            lambda row: deprot_encoded_dictionary.get((row[f'{level}_deprotection']), None), axis=1
        )
    
    # Pull only the encoded columns
    deprot_encoded_columns = select_columns_by_keyword(deprots_96, '_deprot_encoded')

    # Filter the DataFrame to keep only the encoded columns
    deprots_96 = deprots_96[deprot_encoded_columns]

    # Fold the DataFrame to ensure each column has a maximum of 8 values
    # Final DataFrame should have 12 columns, padded if necessary with encoded values
    deprots_96 = fold_dataframe(deprots_96)

    # ================================
    # Functions for Acceptors
    # ================================
    # Inputs:
    # context_history: The context history DataFrame containing acceptor information.   
    # indexed_values: List of index values for tetramers.
    # ================================

    # Create a dimer encoded dictionary
    dimer_encoded_dictionary = {
        'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_SMAC><K63_ABOC>)>)' : 9,
        'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_ABOC><K63_SMAC>)>)' : 10,
        'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_ABOC><K63_ABOC>)>)' : 11,

        'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_SMAC>)><K63_ABOC>)': 12,
        'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_SMAC><K63_ABOC>)><K63_ABOC>)': 13,
        'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_ABOC><K63_SMAC>)><K63_ABOC>)': 14,

        'his-GG-1ubq-1-(<K63_1ubq-2-(<K48_SMAC><K63_ABOC>)>)': 15,
        'his-GG-1ubq-1-(<K63_1ubq-2-(<K48_ABOC><K63_SMAC>)>)': 16,
        'his-GG-1ubq-1-(<K63_1ubq-2-(<K48_ABOC><K63_ABOC>)>)': 17,

        'his-GG-1ubq-1-(<K48_ABOC><K63_1ubq-2-(<K63_SMAC>)>)': 18,
        'his-GG-1ubq-1-(<K48_ABOC><K63_1ubq-2-(<K48_SMAC><K63_ABOC>)>)': 19,
        'his-GG-1ubq-1-(<K48_ABOC><K63_1ubq-2-(<K48_ABOC><K63_SMAC>)>)': 20
    }

    # Filter the combined database for reactions and donors based on indexed values
    # gonna have to be pulled from context_history
    filtered_acceptors = context_history[(context_history['index'].isin(indexed_values))]

    # Preserve order of index_list
    filtered_acceptors_sorted = filtered_acceptors.set_index('index').loc[indexed_values].reset_index()

    # Select columns that contain 'formation' in the name
    dimer_acceptor_columns = select_columns_by_keyword(context_history, 'dimer_formation')
    
    # Get the dimer acceptor formation data
    dimer_acceptors_96 = filtered_acceptors_sorted[dimer_acceptor_columns]

    # Get the multimer string name 
    dimer_acceptors_96['dimer_names'] = dimer_acceptors_96['dimer_formation'].apply(
        lambda x: ast.literal_eval(x).get('multimer_string_name') if isinstance(x, str) else None
        )

    # Create acceptor encoded plate maps
    dimer_acceptors_96['dimers_encoded'] = dimer_acceptors_96['dimer_names'].map(dimer_encoded_dictionary)

    # Filter the DataFrame to keep only the encoded columns
    dimer_acceptors_96 = dimer_acceptors_96[['dimers_encoded']]

    # Fold the DataFrame to ensure each column has a maximum of 8 values
    # Final DataFrame should have 12 columns, padded if necessary with encoded values
    dimer_acceptors_96 = fold_dataframe(dimer_acceptors_96)
   
    # ================================
    # Functions for Padding and Finalization
    # ================================

    output_enzymes_donors_96 = pad_and_validate_dataframe(enzymes_donors_96, required_columns=12)
    output_deprots_96 = pad_and_validate_dataframe(deprots_96, required_columns=12)
    output_dimer_acceptors_96 = pad_and_validate_dataframe(dimer_acceptors_96, required_columns=12)
    
    # ================================
    # Prepare the output dictionary
    # This dictionary will contain the encoded DataFrame and the encoded dictionary
    # It can be used for further processing or analysis.
    # The dimer_acceptors can be used for acceptor plate maps or other analyses.
    # ================================
    output_dict = {
        'enzymes_donors_96': output_enzymes_donors_96,
        'deprots_96': output_deprots_96,
        'dimer_acceptors_96': output_dimer_acceptors_96,
        'e_d_encoded_dictionary': e_d_encoded_dictionary,
        'deprot_encoded_dictionary': deprot_encoded_dictionary,
        'dimer_encoded_dictionary': dimer_encoded_dictionary
    }

    return output_dict
