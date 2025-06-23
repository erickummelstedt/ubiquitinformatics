import json 
import logging
import copy
import sys
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
# Load the configuration, all within the runfile 
# =========================================

# change this to 4 or 5
multimer_size = 4
output_dir = project_root / 'back_end' / 'data' / 'reaction_database' / f'multimer_size_{multimer_size}'
# Read the CSV files into DataFrames
ubiquitin_history = pd.read_csv(output_dir / "ubiquitin_history.csv")
reaction_history = pd.read_csv(output_dir / "reaction_history.csv")
donor_history = pd.read_csv(output_dir / "donor_history.csv")
context_history = pd.read_csv(output_dir / "context_history.csv")


ubiquitin_history, context_history = global_deprotection_dual(ubiquitin_history, context_history)