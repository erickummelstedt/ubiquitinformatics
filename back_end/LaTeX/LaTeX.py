from pathlib import Path
from datetime import date
from pylatex import Document, Section, Subsection, Command
from pylatex.utils import NoEscape

import subprocess

# Compile with pdflatex (must be installed on your system)
subprocess.run(['pdflatex', '/Users/ekummelstedt/le_code_base/ubiquitinformatics/back_end/LaTeX/polyubiquitin_graph_theory.tex'])