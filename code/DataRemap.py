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


def get_target_proteinIDs(hlData, backupPath):
    try:
        with open(backupPath): pass
    except IOError:
        def extract_target_proteins():
            majorIDs, otherData = hlData.get_all_data()
            majorIDs = np.unique(majorIDs[majorIDs != 'None'])
            with open(backupPath, 'w') as f:
                np.savetxt(f, majorIDs, fmt=['%s'])
        extract_target_proteins()

    pIDs = None
    with open(backupPath, 'r') as f:
        pIDs = np.genfromtxt(f, dtype=None)

    return pIDs


# intended to be called by remap() below
def _get_data_of_target(targetProtList, sourceProtIds, sourceData, outputPath=None):
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
                  dtype=[('PrimaryID', targetProtList.dtype.str),
                         ('hasValue',int)])
    structTargetProtList['PrimaryID'] = targetProtList
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


def remap(species):
    #outbase = os.path.join(CLEAN_INPUT_PATH, species.rel_path())
    masterProtListPath = os.path.join(
        CLEAN_INPUT_PATH, 
        species.rel_path('proteins_master_list.txt'))    
    masterProtList = get_target_proteinIDs(
        species.HLDataClass(), masterProtListPath
        )
                    
    for dataCls in species.dataClasses:
        dataObj = dataCls()
        name = re.search(r'(.*)\.\w+', dataObj.featureFileName).group(1)
        outpath = os.path.join(CLEAN_INPUT_PATH,
                               name+'.derived_features.tsv') # name includes path species specific path info
        allIds, allData = dataObj.get_all_data()
        _get_data_of_target(masterProtList, allIds, allData, outpath)


# if __name__ == '__main__':
#     remap(DegRateData())
#     remap(DisorderData()) 
#     remap(SeqData())
#     remap(UbqSitesData())
#     remap(AcetSitesData())
#     remap(MetSitesData())
#     remap(PhosphoSitesData())
#     remap(DegrMotifData())
#     remap(SecStructData())

