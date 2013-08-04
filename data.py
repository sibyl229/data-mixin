from __future__ import division
import os
import numpy as np
import re
from Bio import SeqIO
from abc import ABCMeta, abstractmethod
from local_paths import *

class AbstractData(object):

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


class DisorderData(AbstractData):

    def __init__(self, forceCompute=False):
        self.featureFileName = 'MDB_Homo_sapiens_disorder_consensus.csv'
        self.featureFilePath = os.path.join(RAW_INPUT_PATH, self.featureFileName)
        idColName = 'Uniprot'

        try:
            with open(self.featureFilePath,'r'): pass
        except IOError:
            self.compute_scores()
        else:
            if forceCompute:
                self.compute_scores()

        super(DisorderData, self).__init__(
            self.featureFileName, idColName)


    def compute_scores(self):
        inpath = os.path.join(RAW_INPUT_PATH, 
                              'MDB_Homo_sapiens/disorder_consensus.fasta')

        def score(record):
            aaDisorder = self.parse_disorder_consensus(record.seq) # e.g [1,1,1,0,...1,1], where 1 means disordered
            aaDisorderStr = ''.join(aaDisorder.astype(int).astype(str)) # e.g '1110...11'
            totalDisorderAA = np.sum(aaDisorder)
            ntermDisorder = np.sum(aaDisorder[:40])/40
            internalDisorderRegions = re.findall(r'1+',aaDisorderStr[40:])
          #  import pdb; pdb.set_trace()
            internaDisorderCnt = len([disOR for disOR in internalDisorderRegions if len(disOR)>50])
            return (rec.id, 
                    totalDisorderAA,
                    ntermDisorder,
                    internaDisorderCnt
                )

        dtype=[('Uniprot', '|S20'),
               ('totalDisorderAA', int),
               ('ntermDisorder', float),
               ('internaDisorderCnt', int),
        ]
        with open(inpath, "rU") as handle:
            allScores = [score(rec) for rec in SeqIO.parse(handle, "fasta")]
            allScores = np.array(allScores,dtype=dtype)
        header = [fieldName for fieldName, fieldType in dtype]

        self.backup_scores(allScores,
                           header,
                           self.featureFilePath)

        return allScores, header


    def parse_disorder_consensus(self, seq):
        # get array of amino acid -wise disorder consensus score
        aaDisorder = np.fromiter(seq, dtype=int, count=len(seq))
        return aaDisorder >= 5

        
    def backup_scores(self, scores, header, outpath):
 #       header = [fieldName for fieldName, fieldType in dtype]
 #       fmt = [fieldType for fieldName, fieldType in dtype]
        with open(outpath, 'w') as f:
            delimiter = '\t'
            f.write(delimiter.join(header)+'\n')
            np.savetxt(f, scores, fmt='%s', delimiter=delimiter)
        return True
        
        
if __name__ == '__main__':
    DisorderData(forceCompute=True)
