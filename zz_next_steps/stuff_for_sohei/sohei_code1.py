from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import re
import subprocess

# Load the uploaded Excel file
file_path = Path("/Users/ekummelstedt/Desktop/table_1.xlsx")
df = pd.read_excel(file_path)
# Insert column headers as the first row in the DataFrame
df.loc[-1] = df.columns  # Add header as first row
df.index = df.index + 1  # Shift index
df = df.sort_index()     # Reorder so header is first

# Variables for easy adjustments
var_font_size = 10
var_row_height = 0.8
var_line_width = 0.5

cell_width_x = -0.5
cell_width_y = 0

# Font definitions
font_default = FontProperties(family="DejaVu Sans", size=var_font_size)

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

def format_text(text):
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

# Draw only the middle horizontal line
middle_y = df.shape[0] // 2
ax.plot([-cell_width_y - 0.5, df.shape[1] + cell_width_y + 0.5], [middle_y, middle_y], color='gray', linewidth=0.8)

# Save the new text-based layout
output_path = "/Users/ekummelstedt/Desktop/table_1.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True)

print(f"Saved table image to: {output_path}")

subprocess.run(["open", output_path])  # For macOS