from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import re
import subprocess
import sys

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

# ============================================================

idx = indexes[0] 

def single_dict_for_reaction_schemes(data_dict, idx, multimer_size):
    
    """ 
    Function to create a dictionary for reaction schemes based on the provided index and multimer size.
    This function retrieves the combined database and context history from the data dictionary,
    filters the combined database for the specified index, and extracts relevant columns to form a dictionary.
    
    Parameters:
    
    - data_dict (dict): Dictionary containing the combined database and context history.
    - idx (int): The index for which the reaction scheme is to be created.  
    - multimer_size (int): The size of the multimer (e.g., 4 for tetramers, 5 for pentamers).
    
    Returns:
    
    - steps_dict (list): A list containing a dictionary with the reaction scheme details.
    The dictionary includes the simulation index, multimer ID, acceptor dimer code, and combined values
    for the steps in the reaction scheme.
    The keys in the dictionary are updated to have more readable names for the UI.
    The values are formatted to include the acceptor dimer code and combined values for each step
    in the reaction scheme.
    The function also handles the representation of different steps in the reaction scheme,
    such as deprotection and formation steps, and updates the values accordingly.   
    
    """

    combined_database = data_dict['combined_database']
    context_history = data_dict['context_history']

    # Get the acceptor dimer name from the context history

    def get_acceptor_dimer_name(context_history, idx):
        """ Get the acceptor dimer code from the working DataFrame. """
        
        # Create a dimer encoded dictionary
        dimer_encoded_str_dictionary = {
            # maybe change the code
            'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_SMAC><K63_ABOC>)>)' : "Ub₂ᴬ 9",
            'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_ABOC><K63_SMAC>)>)' : "Ub₂ᴬ 10",
            'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_ABOC><K63_ABOC>)>)' : "Ub₂ᴬ 11",
            'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_SMAC>)><K63_ABOC>)': "Ub₂ᴬ 12",
            'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_SMAC><K63_ABOC>)><K63_ABOC>)': "Ub₂ᴬ 13",
            'his-GG-1ubq-1-(<K48_1ubq-2-(<K48_ABOC><K63_SMAC>)><K63_ABOC>)': "Ub₂ᴬ 14",
            'his-GG-1ubq-1-(<K63_1ubq-2-(<K48_SMAC><K63_ABOC>)>)': "Ub₂ᴬ 15",
            'his-GG-1ubq-1-(<K63_1ubq-2-(<K48_ABOC><K63_SMAC>)>)': "Ub₂ᴬ 16",
            'his-GG-1ubq-1-(<K63_1ubq-2-(<K48_ABOC><K63_ABOC>)>)': "Ub₂ᴬ 17",
            'his-GG-1ubq-1-(<K48_ABOC><K63_1ubq-2-(<K63_SMAC>)>)': "Ub₂ᴬ 18",
            'his-GG-1ubq-1-(<K48_ABOC><K63_1ubq-2-(<K48_SMAC><K63_ABOC>)>)': "Ub₂ᴬ 19",
            'his-GG-1ubq-1-(<K48_ABOC><K63_1ubq-2-(<K48_ABOC><K63_SMAC>)>)': "Ub₂ᴬ 20",
        }

        working_context_df = context_history[context_history['index'] == idx]

        # Get the value of 'context' in 'dimer_formation' column for the selected index
        if not working_context_df.empty:
            context_dict = working_context_df.iloc[0]['dimer_formation']
            # If the value is a string, try to parse it as a dictionary
            if isinstance(context_dict, str):
                import ast
                try:
                    context_dict = ast.literal_eval(context_dict)
                except Exception:
                    context_dict = {}
            multimer_name = context_dict.get('multimer_string_name', None)
            # Map the multimer name to the dimer code using the dictionary
            dimer_code = dimer_encoded_str_dictionary.get(multimer_name, None)
            return dimer_code
        else:
            ValueError(f"No context row found for index: {idx}")
            return None

    # Get the acceptor dimer name for the given index
    multimer_code = get_acceptor_dimer_name(context_history, idx)

    # ============================================================

    # Filter combined_database for rows where index == idx and table_origin is 'Donors' or 'Reactions'
    filtered_combined_df = combined_database[(combined_database['index'] == idx) & (combined_database['table_origin'].isin(['Donors', 'Reactions']))]

    # Determine the last column name based on multimer size
    if multimer_size == 4:
        last_col = 'tetramer_formation'
    elif multimer_size == 5:
        last_col = 'pentamer_formation'
    else:
        raise ValueError(f"Unsupported multimer_size: {multimer_size}")

    # Get all column names
    all_cols = list(filtered_combined_df.columns)

    # Find the start and end indices for slicing columns
    try:
        start_idx = all_cols.index('dimer_deprotection')
        end_idx = all_cols.index(last_col)
    except ValueError as e:
        raise ValueError(f"Column not found: {e}")

    # Slice the columns between dimer_deprotection and last_col (inclusive)
    selected_cols = all_cols[start_idx:end_idx+1]

    # Get the relevant part of the DataFrame
    selected_steps_df = filtered_combined_df[selected_cols]

    # Combine Donors and Reactions rows to form combined_values like "Ub^D 2\ngp78-Ube2g2"
    combined_values = []
    if not selected_steps_df.empty and len(selected_steps_df) == 2:
        # Assume first row is Donors, second is Reactions (order from filter)
        donors_row = selected_steps_df.iloc[0]
        reactions_row = selected_steps_df.iloc[1]
        for donor_val, reaction_val in zip(donors_row, reactions_row):
            combined = f"{donor_val}\n{reaction_val}"
            combined_values.append(combined)
    else:
        ValueError("selected_steps_df does not have exactly 2 rows (Donors and Reactions). Cannot combine.")

    # This is just syntax for the UI
    representation_changes = {
            'Ubc13/Mms2' : 'Ubc13-Mms2',
            'gp78/Ube2g2' : 'gp78-Ube2g2',
            'ubi_ubq_1_K48_ABOC_K63_ABOC' : 'Ubᴰ 5',
            'ubi_ubq_1_K48_SMAC_K63_ABOC' : 'Ubᴰ 2',
            'ubi_ubq_1_K48_ABOC_K63_SMAC' : 'Ubᴰ 4',
            'ubi_ubq_1_K48_SMAC' : 'Ubᴰ 1',
            'ubi_ubq_1_K63_SMAC' : 'Ubᴰ 3',
            'nan\nSMAC_deprot' : 'Smac',
            'nan\nFAKE_deprot' : 'Fake', 
            'dimer_deprotection' : 'Dimer\ndeprotection',
            'trimer_deprotection' : 'Trimer\ndeprotection',
            'tetramer_deprotection' : 'Tetramer\ndeprotection',
            'pentamer_deprotection' : 'Pentamer\ndeprotection',
            'trimer_formation' : 'Trimer\nformation',
            'tetramer_formation' : 'Tetramer\nformation',
            'pentamer_formation' : 'Pentamer\nformation',
        }

    def update_values(values, representation_changes):
        """
        Updates a list of values by replacing exact matches or substrings
        based on a given mapping.

        Parameters:
        - values (list of str): The original values.
        - representation_changes (dict): Mapping of old -> new representations.

        Returns:
        - list of str: Updated values.
        """
        updated_values = []

        for val in values:
            if val in representation_changes:
                updated_val = representation_changes[val]
            else:
                updated_val = val
                for old, new in representation_changes.items():
                    if old in updated_val:
                        updated_val = updated_val.replace(old, new)
            updated_values.append(updated_val)

        return updated_values


    updated_values = update_values(combined_values, representation_changes)
    headers_for_dict = update_values(selected_steps_df.columns, representation_changes)


    # Get the unique multimer_id value from filtered_combined_df
    multimer_id = filtered_combined_df['multimer_id'].unique()[0]

    def format_multimer_id(multimer_id):
        """
        Reformat multimer_id string from 'UbX_Y' to 'Ub_X Y', with suffixes for certain sizes.
        For size 4: add ^T as superscript. For size 5: add ^P as superscript.
        """
        match = re.match(r'^Ub(\d+)_(\d+)$', multimer_id)
        if match:
            size, number = match.group(1), match.group(2)
            suffix = '^T' if size == '4' else '^P' if size == '5' else ''
            return f"Ub_{size}{suffix} {number}"
        return multimer_id
    multimer_id = format_multimer_id(multimer_id)

    # Create a dictionary mapping columns to updated values, with 'Acceptor' as the first key
    steps_dict = {
        'Simulation\nindex': idx, 
        'Multimer Id': multimer_id,
        'Acceptor': multimer_code
        }

    # Update the steps_dict with the combined values
    steps_dict.update(dict(zip(headers_for_dict, updated_values)))
    
    # List the steps_dict to be returned
    steps_dict = [steps_dict]
    # Ensure 'Reaction Number' is always the first key in each dictionary before converting to DataFrame
    steps_dict[0] = {"Reaction\nNumber": None, **steps_dict[0]}

    return steps_dict

def full_dict_df_for_reaction_schemes(data_dict, indexes, multimer_size):
    """
    Function to create a dictionary for reaction schemes based on the provided index and multimer size.
    This function retrieves the combined database and context history from the data dictionary,
    filters the combined_database for the specified index, and extracts relevant columns to form a dictionary.
    """

    steps_full_dict = []
    for i in range(len(indexes)):
        idx = indexes[i]
        # Call the placeholder function to get the steps_dict for each index
        # This is just a placeholder function to avoid syntax errors
        # In practice, this would be replaced with the actual function that processes the data
        steps_dict = single_dict_for_reaction_schemes(data_dict, idx, multimer_size)
        steps_full_dict += steps_dict

    print("steps_full_dict:", steps_full_dict)
    # Add reaction numbers to each step in steps_full_dict
    for i, step in enumerate(steps_full_dict, start=1):
        print("step:", step)
        step['Reaction\nNumber'] = str(i)

    # Create DataFrame preserving original column names
    steps_df = pd.DataFrame(steps_full_dict)

    # Insert headers as the first row
    header_row = {col: col for col in steps_df.columns}
    steps_df = pd.concat([pd.DataFrame([header_row]), steps_df], ignore_index=True)

    return steps_df, steps_full_dict
# ============================================================
# ============================================================

def format_text(text):
    # Unicode superscript/subscript mapping
    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    superscript_map = {
        '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
        '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
        'a': 'ᵃ', 'b': 'ᵇ', 'c': 'ᶜ', 'd': 'ᵈ', 'e': 'ᵉ',
        'f': 'ᶠ', 'g': 'ᵍ', 'h': 'ʰ', 'i': 'ⁱ', 'j': 'ʲ',
        'k': 'ᵏ', 'l': 'ˡ', 'm': 'ᵐ', 'n': 'ⁿ', 'o': 'ᵒ',
        'p': 'ᵖ', 'r': 'ʳ', 's': 'ˢ', 't': 'ᵗ', 'u': 'ᵘ',
        'v': 'ᵛ', 'w': 'ʷ', 'x': 'ˣ', 'y': 'ʸ', 'z': 'ᶻ',
        'A': 'ᴬ', 'B': 'ᴮ', 'D': 'ᴰ', 'E': 'ᴱ', 'G': 'ᴳ',
        'H': 'ᴴ', 'I': 'ᴵ', 'J': 'ᴶ', 'K': 'ᴷ', 'L': 'ᴸ',
        'M': 'ᴹ', 'N': 'ᴺ', 'O': 'ᴼ', 'P': 'ᴾ', 'R': 'ᴿ',
        'T': 'ᵀ', 'U': 'ᵁ', 'V': 'ⱽ', 'W': 'ᵂ'
    }

    # Convert _<digit> to subscript
    text = re.sub(r'_(\d)', lambda m: m.group(1).translate(subscript_map), text)
    # Convert ^<char> to superscript
    text = re.sub(r'\^([A-Za-z0-9])', lambda m: superscript_map.get(m.group(1), m.group(1)), text)
    # Wrap any number immediately following a superscript with {BOLD} marker
    superscript_chars = ''.join(superscript_map.values())
    # Matches a superscript character from your superscript_chars.
	# Ensures it’s followed by a space and digits.
	# Only wraps the digits (\2) in {BOLD}, not the superscript or the space.    
    text = re.sub(r'([' + re.escape(superscript_chars) + r']) (\d+)', r'\1 {BOLD}\2{BOLD}', text)    
    return text

def create_table_image(df):
    """ Create a table image from a DataFrame with custom formatting. """
    
    # Variables for easy adjustments
    var_font_size = 10
    var_row_height = 0.8
    var_line_width = 0.5

    cell_width_x = -0.5
    cell_width_y = 0

    # Font definitions
    font_default = FontProperties(family="DejaVu Sans", size=var_font_size)

    
    # Check if DataFrame is empty
    if df.empty:
        print("DataFrame is empty. No image created.")
        return 
    # Create new plot using plt.text()
    fig, ax = plt.subplots(figsize=(df.shape[1] * 1.5, df.shape[0] * 0.6))
    ax.set_xlim(-cell_width_x - 0.18, df.shape[1] + cell_width_x + 0.4)
    ax.set_ylim(-cell_width_y, df.shape[0] + cell_width_y)
    ax.axis('off')

    # Draw grid-like layout manually with text
    for row in range(df.shape[0]):
        for col in range(df.shape[1]):
            text = str(df.iat[row, col])
            formatted = format_text(text)
            # Check for {BOLD} markers
            bold_parts = re.findall(r'{BOLD}(.*?){BOLD}', formatted)
            clean_text = re.sub(r'{BOLD}(.*?){BOLD}', r'\1', formatted)
            # Add small padding to the text position for clearer separation
            x_offset = col + 0.675 if col == 0 else col + 0.5
            
            # ===================NEW CODE========================
            # Split formatted text into segments: normal and bold
            segments = re.split(r'({BOLD}.*?{BOLD})', formatted)
            has_bold = any(s.startswith("{BOLD}") and s.endswith("{BOLD}") for s in segments)

            if has_bold:
                first_text = re.sub(r'{BOLD}(.*?){BOLD}', r'\1', segments[0])
                if first_text.startswith("Ub₂ᴬ"):
                    cursor_x = x_offset - 0.075
                elif first_text.startswith("Ubᴰ"):
                    cursor_x = x_offset - 0.03
                else:
                    cursor_x = x_offset

                for i, segment in enumerate(segments):
                    if segment.startswith("{BOLD}") and segment.endswith("{BOLD}"):
                        segment_text = segment[6:-6]
                        font = FontProperties(family="DejaVu Sans", size=var_font_size, weight="bold")
                    else:
                        segment_text = segment
                        font = font_default

                    if (len(segments) == 3) & (segments[-1] != ""):
                        if i < 2:
                            y_offset = df.shape[0] - row - 0.38
                        else:
                            y_offset = df.shape[0] - row - 0.62
                    else:
                        y_offset = df.shape[0] - row - 0.5

                    if (i == 2) & (segment != ""):
                        cursor_x = x_offset

                    ax.text(cursor_x, y_offset, segment_text,
                            ha='center', va='center', fontproperties=font)
                    cursor_x += 0.04 * len(segment_text)
            else:
                ax.text(x_offset, df.shape[0] - row - 0.5, clean_text,
                        ha='center', va='center', fontproperties=font_default)
            # ===================       ========================

    # Draw horizontal line below the title row
    title_row_y = df.shape[0] 
    ax.plot([-cell_width_y - 0.5, df.shape[1] + cell_width_y + 0.5], [title_row_y - 1, title_row_y - 1], color='gray', linewidth=0.8)

    # Draw vertical line after the 'Multimer Id' column (if present)
    if 'Multimer Id' in df.columns:
        multimer_col_idx = list(df.columns).index('Multimer Id')
        x_vert = multimer_col_idx + 1  # After the column
        ax.plot([x_vert, x_vert], [0, df.shape[0]], color='gray', linewidth=0.8)

    # Save the new text-based layout as PNG and also get image bytes
    output_path = "/Users/ekummelstedt/Desktop/table_1.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)

    # Save to a BytesIO buffer
    import io
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight', transparent=True)
    buf.seek(0)
    img_bytes = buf.getvalue()
    buf.close()

    # Create a PDF document and embed the image
    from reportlab.lib.pagesizes import letter
    from reportlab.pdfgen import canvas
    from reportlab.lib.utils import ImageReader

    pdf_path = "/Users/ekummelstedt/Desktop/table_1.pdf"
    c = canvas.Canvas(pdf_path, pagesize=letter)
    width, height = letter
    image = ImageReader(io.BytesIO(img_bytes))
    # Fit image to page width, keep aspect ratio
    img_width, img_height = image.getSize()
    aspect = img_height / float(img_width)
    new_width = width - 72  # 1 inch margin on each side
    new_height = new_width * aspect
    c.drawImage(image, 36, height - new_height - 36, width=new_width, height=new_height)
    c.showPage()
    c.save()

    print(f"Saved table image to: {output_path}")
    print(f"Saved PDF to: {pdf_path}")
    return img_bytes, pdf_path


# ============================================================

steps_df, _ = full_dict_df_for_reaction_schemes(data_dict, indexes, multimer_size)

img_bytes, pdf_path = create_table_image(steps_df)

subprocess.run(["open", pdf_path])  # For macOS



# =============================================================

[
  {
    "Unnamed: 0":35,
    "Acceptor":"Ub_2^A 13",
    "step 1":"Smac deprotection",
    "step 2":"Ub^D 2\ngp78-Ube2g2",
    "step 3":"Smac deprotection",
    "step 4":"Ub^D 2\ngp78-Ube2g2",
    "step 5":"Smac deprotection",
    "step 6":"Ub^D 5\ngp78-Ube2g2"
  }
]

# ============================================================