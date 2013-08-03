import numpy as np 
import os
import re
from helpers import ColumnStack
from local_paths import *

'''Take a delimited file of a list of protein information, and select and reorder it,
so that they all have 
    1) the same proteins in the same ordering
    2) without any duplicate''' 


class ReMap(object):

    @classmethod
    def get_data_of_target(cls, protIds, dataArray, outputPath=None):
        targetProtList = cls.get_target_proteinIDs()

        pIDtoIndex = {pid: i \
            for (i,pid) in reversed(list(enumerate(protIds)))} # revered iteration to ensure using first entry, in case of duplicate protein ID among entries
       
        # add an empty tuple to the end of the dataArray, for protein without data
        emptyEntry = tuple(np.zeros(len(dataArray[0])))
        newEndIndex = len(dataArray)
        np.resize(dataArray, newEndIndex+1)
        dataArray[-1] = emptyEntry # I intend to modify the dataArray itself, the warning is puzzling

        @np.vectorize
        def locateEntry(id):
            entryIndex = pIDtoIndex.get(id, newEndIndex)
            return entryIndex

        selectedRows = locateEntry(targetProtList)
        targetProtData = dataArray[selectedRows]

        structTargetProtList = np.zeros((targetProtList.shape[0],),dtype=[('Uniprot', targetProtList.dtype.str)] )
        structTargetProtList['Uniprot'] = targetProtList
        targetIDAndData = ColumnStack.stack(structTargetProtList, targetProtData)
        if outputPath:
            with open(outputPath, 'w') as f:
                header = '\t'.join(targetIDAndData.dtype.names)+'\n'
                f.write(header)
                np.savetxt(f, targetIDAndData, fmt='%s')
     
        return targetIDAndData


    @staticmethod
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
                    


class DegRateData(object):
    inputName = 'pr101183k_si_002_HeLa.csv'
    rawDataPath = os.path.join(RAW_INPUT_PATH, inputName)
    rawdata = np.genfromtxt(rawDataPath, delimiter='\t', 
        usecols=[3,10], 
        dtype=None, names=True)
    idColName = 'Uniprot'
    dataColNames = list(set(rawdata.dtype.names) - {idColName})

    _protIds, _dataArray = None, None

    @classmethod
    def get_data_of_target(cls):
        protIds, dataArray = cls.get_all_data()
        outpath = os.path.join(FEATURE_FILE_PATH, 'featureFrom_'+cls.inputName)
        return ReMap.get_data_of_target(protIds, dataArray, outpath)

    @classmethod
    def get_all_data(cls):
        if cls._protIds == None or cls._dataArray == None:
            rawProtIDs = cls.rawdata[cls.idColName]
            vectorized_get_major_id = np.vectorize(cls.get_major_id)
            cls._protIds = vectorized_get_major_id(rawProtIDs)
            cls._dataArray = cls.rawdata[cls.dataColNames]
        return cls._protIds, cls._dataArray


    @staticmethod
    def get_major_id(idstring):
        match = re.search(r'(\w+)',idstring)
        if match:
            return match.group(1)
        else:
            return None

DegRateData.get_data_of_target()



        # includeColumn='all', delimiter='\t', preprocess=None):
        # if preprocess:
        #     f = open(file_path, 'r')
        #     data = preprocess(file_path)
        # else:
        #     np.genfromtxt()
