from __future__ import division
import os
import numpy as np
import re
from abc import ABCMeta, abstractmethod
from settings import *

class AbstractData(object):
    '''Takes care of reading input raw data files, 
and Used by Remap(), which calls get_all_data().
Intended to be inherited by concrete data classes that deal with specific data. For an example, refer to DegRateData class below'''

    __metaclass__ = ABCMeta
    
    @abstractmethod
    def __init__(self, featureFileName, idColName, inpath=None, usecols=None):
        # read data from file
        self.featureFileName = featureFileName
        inpath = inpath or RAW_INPUT_PATH
        self.featureFilePath = os.path.join(inpath, self.featureFileName)
        self.rawdata = np.genfromtxt(self.featureFilePath, delimiter='\t', 
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
        '''choose first id if multiple exists'''
        match = re.search(r'(\w+)',idstring)
        if match:
            return match.group(1)
        else:
            return None


class DegRateData(AbstractData):

    def __init__(self):
        featureFileName = 'pr101183k_si_002_HeLa.csv'
        idColName = 'Uniprot'
        usecols = [3,11]

        self.featureFileName = featureFileName
        super(DegRateData, self).__init__(
            featureFileName, idColName, usecols=usecols)

