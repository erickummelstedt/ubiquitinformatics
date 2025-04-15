import json 
import logging
import copy
from typing import Dict, List, Any, Union
import sys

# figure out the path issues
# home_dir = os.path.expanduser('~')
# local_path = '/home/erickummelstedt/lecodebase/ubiquitinformatics/src/main.py'
local_path = '/Users/ekummelstedt/le_code_base/ubiquitinformatics'
sys.path.insert(0, local_path)

import pytest
import copy

from src.main import \
    find_branching_site, \
    iterate_through_ubiquitin

