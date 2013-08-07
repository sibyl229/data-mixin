from __future__ import division
import os
import numpy as np
import re
from Bio import SeqIO
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

class AbstractComputedData(object):
    '''Similar to AbstractData class, AbstractComputedData also implements get_all_data() used by remap().
Different from AbstractData class, this one is extended by classes that compute feature scores 
(and backs them up in a feature file).
 '''
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def __init__(self, inputFileName, featureFileName, inpath=None, forceCompute=False):
        inpath = inpath or RAW_INPUT_PATH
        self.rawInputFilePath = os.path.join(inpath, inputFileName)
        
        self.featureFileName = featureFileName
        self.featureFilePath = os.path.join(RAW_INPUT_PATH, self.featureFileName)
 
        try:
            with open(self.featureFilePath,'r'): pass
        except IOError:
            self.compute_and_backup_scores()
        else:
            if forceCompute:
                self.compute_and_backup_scores()

        features = np.genfromtxt(self.featureFilePath, delimiter='\t', 
                dtype=None, names=True)
        self._protIds = features['Uniprot']
        self._dataArray = features[[f for f in features.dtype.names if f != 'Uniprot']]


    def get_all_data(self):
        return self._protIds, self._dataArray

    def compute_and_backup_scores(self):
        scores = self.compute_scores()
        self.backup_scores(scores, self.featureFilePath)
        return scores

    @abstractmethod
    def compute_scores(self):
        '''return a structured numpy array (including field/column names) that includes:
        a 'Uniprot' field for protein id, and additional fields for feature scores'''
        pass
        
    @staticmethod
    def backup_scores(scores, backupPath):
        header = scores.dtype.names
        with open(backupPath, 'w') as f:
            delimiter = '\t'
            f.write(delimiter.join(header)+'\n')
            np.savetxt(f, scores, fmt='%s', delimiter=delimiter)
        return True


class DisorderData(AbstractComputedData):

    def __init__(self, forceCompute=False): 
        rawInputName = 'MDB_Homo_sapiens/disorder_consensus.fasta'
        featureFileName = 'MDB_Homo_sapiens_disorder_consensus.csv'
        super(DisorderData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)


    def compute_scores(self):

        def score(record):
            aaDisorder = self.parse_disorder_consensus(record.seq) # e.g [1,1,1,0,...1,1], where 1 means disordered
            aaDisorderStr = ''.join(aaDisorder.astype(int).astype(str)) # e.g '1110...11'
            totalDisorderAA = np.sum(aaDisorder)
            totalDisorderRatio = totalDisorderAA / len(record.seq)
            ntermDisorder = np.sum(aaDisorder[:40])/40
            internalDisorderRegions = re.findall(r'1+',aaDisorderStr[40:])
          #  import pdb; pdb.set_trace()
            internaDisorderCnt = len([disOR for disOR in internalDisorderRegions if len(disOR)>50])

            return (record.id, 
                    totalDisorderAA,
                    totalDisorderRatio,
                    ntermDisorder,
                    internaDisorderCnt
                )

        dtype=[('Uniprot', '|S20'),
               ('totalDisorderAA', int),
               ('totalDisorderRatio', float),
               ('ntermDisorder', float),
               ('internaDisorderCnt', int),
        ]

        with open(self.rawInputFilePath, "rU") as handle:
#            import pdb; pdb.set_trace()
            allScores = [score(rec) for rec in SeqIO.parse(handle, "fasta")]
            allScores = np.array(allScores,dtype=dtype)

        return allScores


    def parse_disorder_consensus(self, seq):
        # get array of amino acid -wise disorder consensus score
        aaDisorder = np.fromiter(seq, dtype=int, count=len(seq))
        return aaDisorder >= 5
        
        
if __name__ == '__main__':
    DisorderData(forceCompute=True)
