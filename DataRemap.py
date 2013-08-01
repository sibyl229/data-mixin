import numpy as np 
import os
import re

'''Take a delimited file of a list of protein information, and select and reorder it,
so that they all have 
    1) the same proteins in the same ordering
    2) without any duplicate''' 

CURR_PATH = os.path.dirname(__file__)
RAW_INPUT_PATH = os.path.join(CURR_PATH, 'raw_inputs/')

class ReMap(object):

    @classmethod
    def get_data_of_target(cls, protIds, dataArray):
        targetProtList = cls.get_target_proteinIDs()

        pIDtoIndex = {pid: i \
            for (i,pid) in reversed(list(enumerate(protIds)))} # revered iteration to ensure using first entry, in case of duplicate protein ID among entries
        
        def locateEntry(id):
            entryIndex = pIDtoIndex.get(id, -1)
            if entryIndex == -1:
                entry = ['' for c in otherColNums]
            else:
                entry = list(dataArray[entryIndex])
                
            import pdb; pdb.set_trace()
            return entry
        vectorizedLocateEntry = np.vectorize(locateEntry,otypes=[object])

        selectedRowData = locateEntry(targetProtList)
        return selectedRowData


    @staticmethod
    def get_target_proteinIDs():    
        proteinsPath = os.path.join(CURR_PATH, 'input/', 'target_proteins.txt')
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

    rawDataPath = os.path.join(RAW_INPUT_PATH, 'pr101183k_si_002_HeLa.csv')
    rawdata = np.genfromtxt(rawDataPath, delimiter='\t', 
        usecols=[3,11], 
        dtype=None, names=True)
    idColName = 'Uniprot'
    dataColNames = list(set(rawdata.dtype.names) - {idColName})

    _protIds, _dataArray = None, None

    @classmethod
    def get_data_of_target(cls):
        protIds, dataArray = cls.get_all_data()
        return ReMap.get_data_of_target(protIds, dataArray)

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
        print idstring
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
