import numpy as np
import os
import json
from settings import *

'''deals with the use of GFF (general feature format).
Handles retrival of corresponding rows based on proteins'''

gffName = 'uniprot-taxonomy%3A9606+AND+reviewed%3Ayes+AND+keyword%3A181.gff'
gffPath = os.path.join(RAW_INPUT_PATH, gffName)
gffHeader = ['Uniprot','source','feature','start','end','attr']
locatorPath = os.path.join(RAW_INPUT_PATH, 
                           gffName+'.json')

def setup_gf():
    gff = np.genfromtxt(gffPath, dtype=None, 
                  usecols=[0,1,2,3,4,8],
                  missing='.', delimiter='\t',
                  names=gffHeader)

    try:
        with open(locatorPath, 'r') as f:
            locate = json.load(f)
    except IOError:
        locate = {}
        for i,entry in enumerate(gff):
            pid = entry[0]
            if not locate.get(pid):
                locate[pid] = []
            locate[pid].append(i)
        with open(locatorPath, 'w') as f:
            f.write(json.dumps(locate))

    return gff, locate

gff, gff_locate = setup_gf()


def get_gf(protId):
    indices = gff_locate.get(protId, [])
    return gff[indices]



def get_longest_chain(protId):
    '''get information about the longest Chain of the mature protein'''
    protGf = get_gf(protId)
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
