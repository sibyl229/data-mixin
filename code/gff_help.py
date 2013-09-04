import numpy as np
import os
import json
from settings import *

'''deals with the use of GFF (general feature format).
Handles retrival of corresponding rows based on proteins'''

class Gff(object):

    gffHeader = ['Uniprot','source','feature','start','end','attr']
    
    def __init__(self, species):
        gffName = species.context['gffName']
        self.gffPath = os.path.join(RAW_INPUT_PATH, 
                                    species.rel_path(gffName))
        self.locatorPath = os.path.join(RAW_INPUT_PATH, 
                                        species.rel_path(gffName+'.json'))
        self.setup_gf()

    def setup_gf(self):
        gff = np.genfromtxt(self.gffPath, dtype=None, 
                      usecols=[0,1,2,3,4,8],
                      missing_values='.', delimiter='\t',
                      names=self.gffHeader)

        try:
            with open(self.locatorPath, 'r') as f:
                locate = json.load(f)
        except IOError:
            self.locate = {}
            for i,entry in enumerate(gff):
                pid = entry[0]
                if not locate.get(pid):
                    locate[pid] = []
                locate[pid].append(i)
            with open(locatorPath, 'w') as f:
                f.write(json.dumps(locate))

        self.gff = gff
        self.gff_locate = locate
        return


    def get_gf(self, protId):
        indices = self.gff_locate.get(protId, [])
        return self.gff[indices]


    def get_longest_chain(self, protId):
        '''get information about the longest Chain of the mature protein'''
        protGf = self.get_gf(protId)
        isChain = protGf['feature']=='Chain'
        protGf = protGf[isChain]
        if len(protGf) > 0:
            chainLengths = protGf['end'] - protGf['start']
            indexLongest = np.argmax(chainLengths)
            start, end = protGf[['start', 'end']][indexLongest]

            # change 1-indexed position to 0-indexed position
            # AND change closed range [start, end] to open range [start, end)
            # so they follow programming conventions
            start, end = start-1, end
            return start, end
        else:
            return None

# import pdb; pdb.set_trace() 
