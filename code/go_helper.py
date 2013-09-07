import numpy as np
import os
from settings import *

class GO(object):    

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
    
    def __init__(self, species):
        goaName = species.context['goaName']
        self.taxonCode = species.context.get('taxonNum', None)
        self.goPath = os.path.join(RAW_INPUT_PATH, 
                                    species.rel_path(goaName))
        self.cytoplasmicProteins = self.find_cytoplasmic()

    
    def find_cytoplasmic(self):
        goData = np.genfromtxt(self.goPath, 
                                    dtype=None, 
                                    delimiter='\t', 
                                    skip_header=31,
                                    comments='!',
                                    names=self.goDataHeader)
        
        relevant = np.logical_and(goData['taxon'] == 'taxon:%s' % self.taxonCode,
                                  goData['evidence'] != 'IEA')
        inCytoplasm = np.logical_and(
            relevant,
            np.logical_or(goData['GOID'] == 'GO:0005829',
                          goData['GOID'] == 'GO:0005737')
        )
        cytoplasmic = goData[inCytoplasm]
        protIDs = set(cytoplasmic['objID'])
        protIDs.update(cytoplasmic['objSynonym']) 
        #so multiple ids of the protein can be used to determine if its cytoplasmic
        return protIDs

    def is_cytoplasmic(self, pID):
        return pID in self.cytoplasmicProteins



