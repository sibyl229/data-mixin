'''path to several useful folders'''
import os
from helpers import make_sure_path_exists

PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../'))
RAW_INPUT_PATH = os.path.join(PROJECT_ROOT, 'raw_inputs/')
CLEAN_INPUT_PATH = os.path.join(PROJECT_ROOT, 'input/')
FEATURE_FILE_PATH = os.path.join(PROJECT_ROOT, 'input/features')
#NORMALIZED_FEATURE_FILE_PATH = os.path.join(PROJECT_ROOT, 'input/normalized_features')
FIG_PATH = os.path.join(PROJECT_ROOT, 'fig')

for p in [RAW_INPUT_PATH, CLEAN_INPUT_PATH, FIG_PATH]:
    make_sure_path_exists(p)
