from pathlib import Path
import pandas as pd
import sys

# Dynamically get the backend path relative to this file
current_file = Path(__file__).resolve()
project_root = current_file.parents[2]  # Go up to project root
sys.path.insert(0, str(project_root))
local_path = project_root / 'back_end'
sys.path.insert(0, str(local_path))

from src.simulation import create_reaction_histories
from src.utils.utils import get_multimer_column_names
from tests.test_data import (
    histag_ubi_ubq_1,
    histag_ubi_ubq_1_K48_aboc,
    histag_ubi_ubq_1_K63_aboc,
    ubi_ubq_1_K48_SMAC,
    ubi_ubq_1_K63_SMAC,
    ubi_ubq_1_K48_SMAC_K63_ABOC,
    ubi_ubq_1_K48_ABOC_K63_SMAC,
    ubi_ubq_1_K48_ABOC_K63_ABOC
)

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
for i in range(2,7): 
    save_reaction_database(
        acceptor_list=acceptor_list, 
        donor_list=donor_list, 
        multimer_size=i
    )
