'''path to several useful folders'''
import os

PROJECT_ROOT = os.path.dirname(__file__)
RAW_INPUT_PATH = os.path.join(PROJECT_ROOT, 'raw_inputs/')
CLEAN_INPUT_PATH = os.path.join(PROJECT_ROOT, 'input/')
FEATURE_FILE_PATH = os.path.join(PROJECT_ROOT, 'input/features')
NORMALIZED_FEATURE_FILE_PATH = os.path.join(PROJECT_ROOT, 'input/normalized_features')
