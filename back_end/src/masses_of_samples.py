import ast
import pandas as pd
import sys
from pathlib import Path

# Dynamically get the backend path relative to this file
current_file = Path(__file__).resolve()
project_root = current_file.parents[2]  # Go up to project root
sys.path.insert(0, str(project_root))
local_path = project_root / 'back_end'
sys.path.insert(0, str(local_path))

from src.main import iterate_through_ubiquitin
from src.utils.utils import convert_json_to_dict


# Protecting group masses
MASSES = {
    'ABOC': 131.13,  # Allyloxycarbonyl group
    'SMAC': 144.13,  # S-methyl allyloxycarbonyl group
    'UBIQUITIN_DHHHHHH': 9526.81,  # Monoisotopic mass of ubiquitin with His-tag (DHHHHHH)
    'WATER': 18.01056,  # Mass of water (H2O) lost during isopeptide bond formation
    'UBIQUITIN': 8588.87,  # Monoisotopic mass of ubiquitin with acetylated N-terminus with M to Nle mutation (no His-tag)
}

# ------------------------------------------------------------------
# Data loading
# ------------------------------------------------------------------

DATA_DIR = project_root / 'back_end' / 'data' / 'filtered_reaction_database'

UBIQUITIN_HISTORY_4 = DATA_DIR / 'multimer_size_4' / 'ubiquitin_history.csv'
UBIQUITIN_HISTORY_5 = DATA_DIR / 'multimer_size_5' / 'ubiquitin_history.csv'


def load_ubiquitin_history(csv_path):
    """Load a ubiquitin_history CSV into a DataFrame, filtered to used_in_synthesis == 1."""
    df = pd.read_csv(csv_path)
    df = df[df['used_in_synthesis'] == 1].reset_index(drop=True)
    return df


def parse_ubiquitin_dict(cell_value):
    """Safely parse a stringified Python dict from a CSV cell."""
    if pd.isna(cell_value) or cell_value == '':
        return None
    return ast.literal_eval(cell_value)


def process_ubiquitin_history(df):
    """
    For each row and each synthesis step column, parse the ubiquitin dict
    and run iterate_through_ubiquitin to get the relabeled dict and context.

    Returns a dict of {(row_index, column_name): (adapted_dict, context)}
    """
    TETRAMER_STEPS = [
        'initial_acceptor', 'dimer_formation', 'dimer_deprotection',
        'trimer_formation', 'trimer_deprotection', 'tetramer_formation',
    ]
    PENTAMER_STEPS = TETRAMER_STEPS + ['tetramer_deprotection', 'pentamer_formation']

    # Use whichever steps are present in the DataFrame
    step_columns = [c for c in PENTAMER_STEPS if c in df.columns]

    results = {}
    for row_idx, row in df.iterrows():
        for col in step_columns:
            ub_dict = parse_ubiquitin_dict(row[col])
            if ub_dict is None:
                results[(row_idx, col)] = (None, None)
                continue
            adapted_dict, context = iterate_through_ubiquitin(ub_dict)
            results[(row_idx, col)] = (adapted_dict, context)
    return results


def calculate_mass(context):
    """
    Calculate the total mass of a ubiquitin multimer from its context.

    Mass = 1 * UBIQUITIN_DHHHHHH + (max_chain_number - 1) * (UBIQUITIN - WATER)
         + len(ABOC_lysines) * ABOC
         + len(SMAC_lysines) * SMAC
    """
    chain_mass = MASSES['UBIQUITIN_DHHHHHH'] + (context['max_chain_number'] - 1) * (MASSES['UBIQUITIN'] - MASSES['WATER'])
    aboc_mass = len(context['ABOC_lysines']) * MASSES['ABOC']
    smac_mass = len(context['SMAC_lysines']) * MASSES['SMAC']
    return chain_mass + aboc_mass + smac_mass


def build_mass_table(df):
    """
    Build a DataFrame of masses for each multimer at each synthesis step.
    """
    results = process_ubiquitin_history(df)

    TETRAMER_STEPS = [
        'initial_acceptor', 'dimer_formation', 'dimer_deprotection',
        'trimer_formation', 'trimer_deprotection', 'tetramer_formation',
    ]
    PENTAMER_STEPS = TETRAMER_STEPS + ['tetramer_deprotection', 'pentamer_formation']
    step_columns = [c for c in PENTAMER_STEPS if c in df.columns]

    rows = []
    for row_idx, row in df.iterrows():
        entry = {'multimer_id': row['multimer_id']}
        for col in step_columns:
            adapted_dict, context = results[(row_idx, col)]
            if context is not None:
                entry[col] = round(calculate_mass(context), 2)
            else:
                entry[col] = None
        rows.append(entry)
    result = pd.DataFrame(rows)
    result['_sort_key'] = result['multimer_id'].str.extract(r'(\d+)$').astype(int)
    result = result.sort_values('_sort_key').drop(columns='_sort_key').reset_index(drop=True)
    return result


if __name__ == '__main__':
    output_path = project_root / 'back_end' / 'data' / 'masses_of_samples.xlsx'

    # Tetramers
    print("Processing tetramers...")
    df4 = load_ubiquitin_history(UBIQUITIN_HISTORY_4)
    mass_table_4 = build_mass_table(df4)
    print(mass_table_4.to_string(index=False))

    # Pentamers
    print("\nProcessing pentamers...")
    df5 = load_ubiquitin_history(UBIQUITIN_HISTORY_5)
    mass_table_5 = build_mass_table(df5)
    print(mass_table_5.to_string(index=False))

    # Write to Excel with separate sheets
    with pd.ExcelWriter(output_path) as writer:
        mass_table_4.to_excel(writer, sheet_name='tetramers', index=False)
        mass_table_5.to_excel(writer, sheet_name='pentamers', index=False)

    print(f"\nExcel file saved to: {output_path}")



