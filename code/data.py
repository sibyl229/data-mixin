from __future__ import division
import os
import numpy as np
import re
from abc import abstractmethod
from settings import *

class AbstractData(object):
    '''Takes care of reading input raw data files, 
and Used by Remap(), which calls get_all_data().
Intended to be inherited by concrete data classes that deal with specific data. For an example, refer to DegRateData class below'''

    
    @abstractmethod
    def __init__(self, featureFileName, idColName, usecols=None, skip_header=0, names=True):
        # read data from file
        self.featureFileName = featureFileName
        self.featureFilePath = os.path.join(RAW_INPUT_PATH, 
                                            self.featureFileName)
        #import pdb; pdb.set_trace()
        self.rawdata = np.genfromtxt(self.featureFilePath, 
                                     delimiter='\t', 
                                     skip_header=skip_header, 
                                     usecols=usecols,  
                                     # Note None would include all columns
                                     dtype=None, 
                                     names=names)

        # extract protein id and features
        self.idColName = idColName
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
        return self._protIds, self._dataArray

    @staticmethod
    def get_major_id(idstring):
        '''choose first id if multiple exists'''
        match = re.search(r'(\w+)',idstring)
        if match:
            return match.group(1)
        else:
            return None

'''
This class is kept as an example. For practical use, include such in submodules of a under package species
class DegRateData(AbstractData):

    def __init__(self):
        featureFileName = 'pr101183k_si_002_HeLa.csv'
        idColName = 'Uniprot'
        usecols = [3,11]

        self.featureFileName = featureFileName
        super(DegRateData, self).__init__(
            featureFileName, idColName, usecols=usecols)

'''
