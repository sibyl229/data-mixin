import numpy as np 
import os
import re
from helpers import *
from settings import *
from data import * 
from compute_data import * # DegRateData

'''Take a delimited file of a list of protein information, and select and reorder it,
so that they all have 
    1) the same proteins in the same ordering
    2) without any duplicate''' 


def get_target_proteinIDs():
    proteinsPath = os.path.join(CLEAN_INPUT_PATH, 'target_proteins.txt')
    try:
        with open(proteinsPath): pass
    except IOError:
        def extract_target_proteins():
            majorIDs, otherData = DegRateData.get_all_data()
            majorIDs = np.unique(majorIDs[majorIDs != 'None'])
            with open(proteinsPath, 'w') as f:
                np.savetxt(f, majorIDs, fmt=['%s'])
        extract_target_proteins()

    pIDs = None
    with open(proteinsPath, 'r') as f:
        pIDs = np.genfromtxt(f, dtype=None)

    return pIDs

targetProtList = get_target_proteinIDs()

# intended to be called by remap() below
def _get_data_of_target(sourceProtIds, sourceData, outputPath=None):
#    targetProtList = get_target_proteinIDs()
    pIDtoIndex = {pid: i \
        for (i,pid) in reversed(list(enumerate(sourceProtIds)))} # revered iteration to ensure using first entry, in case of duplicate protein ID among entries

    # add an empty tuple to the end of the sourceData, for protein without data
    emptyEntry = tuple(np.zeros(len(sourceData[0])))
    newEndIndex = len(sourceData)
    sourceData = np.resize(sourceData, newEndIndex+1)

    sourceData[-1] = emptyEntry # I intend to modify the sourceData itself, the warning is puzzling

    @np.vectorize
    def locateEntry(id):
        entryIndex = pIDtoIndex.get(id, newEndIndex)
        return entryIndex

    selectedRows = locateEntry(targetProtList)
    targetProtData = sourceData[selectedRows]

    structTargetProtList = \
        np.zeros((targetProtList.shape[0],),
                  dtype=[('Uniprot', targetProtList.dtype.str),
                         ('hasValue',int)])
    structTargetProtList['Uniprot'] = targetProtList
    structTargetProtList['hasValue'] = (selectedRows!=newEndIndex)
    targetIDAndData = column_stack(
        structTargetProtList, 
        targetProtData,
    )
    if outputPath:
        with open(outputPath, 'w') as f:
            header = '\t'.join(targetIDAndData.dtype.names)+'\n'
            f.write(header)
            np.savetxt(f, targetIDAndData, fmt='%s', delimiter='\t')

    return targetIDAndData


def remap(dataObj):
    name = re.search(r'(.*)\.\w+', dataObj.featureFileName).group(1)
    outpath = os.path.join(FEATURE_FILE_PATH, 'featureFrom_'+name+'.csv')
    allIds, allData = dataObj.get_all_data()
    return _get_data_of_target(allIds, allData, outpath)


if __name__ == '__main__':
    remap(DegRateData())
    remap(DisorderData()) 
    remap(SeqData())
    remap(UbqSitesData())
    remap(AcetSitesData())
    remap(MetSitesData())
    remap(PhosphoSitesData())


