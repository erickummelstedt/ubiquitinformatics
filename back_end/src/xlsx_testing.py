import json 
import logging
import copy
import sys
import ast
import numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager
import ast 
import subprocess
from collections import defaultdict

# Dynamically get the backend path relative to this script location
script_path = Path(__file__).resolve()
project_root = script_path.parents[2]  # Go up to project root (adjust if needed)
src_path = project_root / 'back_end' / 'src'
sys.path.insert(0, str(src_path))

from utils.utils import *
from utils.logging_utils import *
from main import *
from plotting import *

multimer_size = 4

# ============================================================
# Functions that exist in Fast API 
# ===========================================================

# download CSV files
def download_data_dict(multimer_size):
    input_dir = project_root / 'back_end' / 'data' / 'filtered_reaction_database' / f'multimer_size_{multimer_size}'

    combined_database = pd.read_csv(input_dir / 'combined_database.csv', index_col=0)
    context_history = pd.read_csv(input_dir / 'context_history.csv', index_col=0)
    donor_history = pd.read_csv(input_dir / 'donor_history.csv', index_col=0)
    reaction_history = pd.read_csv(input_dir / 'reaction_history.csv', index_col=0)
    ubiquitin_history = pd.read_csv(input_dir / 'ubiquitin_history.csv', index_col=0)

    return {
        'combined_database': combined_database,
        'context_history': context_history,
        'donor_history': donor_history,
        'reaction_history': reaction_history,
        'ubiquitin_history': ubiquitin_history
    }

# Create plate dataframes for the selected indexes
data_dict = download_data_dict(multimer_size)
combined_database = data_dict['combined_database']
context_history = data_dict['context_history']
donor_history = data_dict['donor_history']
reaction_history = data_dict['reaction_history']
ubiquitin_history = data_dict['ubiquitin_history']

# Create the plate dataframes for the selected indexes
# Select the indexes of the multimer_ids that are used in synthesis
selected_ids = ["Ub4_5","Ub4_10","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_2","Ub4_8", "Ub4_11","Ub4_7"]
indexes = list()

# Get the indexes of the selected multimer_ids that are used in synthesis
for id in selected_ids:
    new_index = int(combined_database[(combined_database['multimer_id'] == id) & (combined_database['used_in_synthesis'] == 1)]['index'].unique()[0])
    indexes.append(new_index)

# Create the plate dataframes for the selected indexes
output_dict = inner_create_plate_dfs(data_dict, indexes, multimer_size)

# start of create_xlsx_bytes function that creates an Excel file with the plate dataframes
excel_bytes = create_xlsx_bytes(output_dict)

# 'excel_bytes' is your BytesIO object
with open("enzyme_donor_and_component_counts.xlsx", "wb") as f:
    f.write(excel_bytes.getbuffer())

# ============================================================
# Close Excel file and open it
# ============================================================

def close_excel_file(filename):
    """
    Ensure the Excel file is closed before proceeding.
    """
    import os
    import time
    # Only attempt on macOS
    if sys.platform != 'darwin':
        return
    # Check if Excel is running
    def is_excel_running():
        try:
            result = subprocess.run([
                'osascript',
                '-e', 'tell application "System Events" to (name of processes) contains "Microsoft Excel"'
            ], capture_output=True, text=True)
            return 'true' in result.stdout.lower()
        except Exception:
            return False
    # Try to close the file if already open (macOS only, using AppleScript)
    if is_excel_running():
        closed = False
        for _ in range(10):  # Try up to 10 times
            try:
                subprocess.run([
                    'osascript',
                    '-e', f'tell application "Microsoft Excel" to close (every document whose name is "{filename}")',
                    '-e', 'delay 0.5'
                ], check=True)
                closed = True
                break
            except Exception:
                time.sleep(0.5)
        if not closed:
            print(f"Warning: Could not close Excel file {filename} automatically.")
        # Wait a moment to ensure file is closed
        time.sleep(1.5)

def close_excel_completely():
    """
    Ensure Microsoft Excel is fully quit before opening the file again.
    """
    import time
    import subprocess
    if sys.platform != 'darwin':
        return
    try:
        subprocess.run([
            'osascript',
            '-e', 'tell application "Microsoft Excel" to quit saving no'
        ], check=True)
        time.sleep(2)  # Wait to ensure Excel has fully quit
    except Exception:
        print("Warning: Could not quit Microsoft Excel automatically.")

close_excel_file('enzyme_donor_and_component_counts.xlsx')
close_excel_completely()
subprocess.run(['open', 'enzyme_donor_and_component_counts.xlsx'])
