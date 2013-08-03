import os
import numpy as np
import re
from abc import ABCMeta, abstractmethod
from local_paths import *

class AbstractData(object):

    __metaclass__ = ABCMeta
    
    @abstractmethod
    def __init__(self, inputName, idColName, inpath=None, usecols=None):
        # read data from file
        self.inputName = inputName
        inpath = inpath or RAW_INPUT_PATH
        rawDataPath = os.path.join(inpath, self.inputName)
        self.rawdata = np.genfromtxt(rawDataPath, delimiter='\t', 
            usecols=usecols, #Note None would include all columns
            dtype=None, names=True)

        # extract protein id and features
        self.idColName = 'Uniprot'
        self.dataColNames = \
            [colName for colName in self.rawdata.dtype.names \
             if colName != self.idColName]
        self._extract_data()
        




    def get_all_data(self):
        return self._protIds, self._dataArray


    def _extract_data(self):
        rawProtIDs = self.rawdata[self.idColName]
        vectorized_get_major_id = np.vectorize(self.get_major_id)
        self._protIds = vectorized_get_major_id(rawProtIDs)
        self._dataArray = self.rawdata[self.dataColNames]
        return

    @staticmethod
    def get_major_id(idstring):
        match = re.search(r'(\w+)',idstring)
        if match:
            return match.group(1)
        else:
            return None

class DegRateData(AbstractData):

    def __init__(self):
        inputName = 'pr101183k_si_002_HeLa.csv'
        idColName = 'Uniprot'
        usecols = [3,11]

        self.inputName = inputName
        super(DegRateData,self).__init__(
            inputName, idColName, usecols=usecols)


