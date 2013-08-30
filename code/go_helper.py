import numpy as np
import os
from settings import *

goPath = os.path.join(RAW_INPUT_PATH,
                      'gene_association.goa_human')

goDataHeader = [
    'DB',
    'objID',
    'objSym',
    'qualifier',
    'GOID',
    'ref',
    'evidence',
    'from',
    'aspect',
    'objName',
    'objSynonym',
    'objType',
    'taxon',
    'date',
    'assignedBy',
    'annoExtension',
    'fromID'
]

goData = np.genfromtxt(goPath, dtype=None, 
                       delimiter='\t', 
                       skip_header=31,
                       comments='!',
                       names=goDataHeader)


relevant = np.logical_and(goData['taxon'] == 'taxon:9606', # human
                          goData['evidence'] != 'IEA')
inCytoplasm = np.logical_and(
    relevant,
    np.logical_or(goData['GOID'] == 'GO:0005829',
                  goData['GOID'] == 'GO:0005737')
)
cytoplasmic = goData[inCytoplasm]
cytoplasmicProteins = set(cytoplasmic['objID'])

@np.vectorize
def is_cytoplasmic(pID):
    return pID in cytoplasmicProteins


