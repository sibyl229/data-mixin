from __future__ import division
import os
import numpy as np
import re
from Bio import SeqIO
from abc import ABCMeta, abstractmethod
from settings import *
from helpers import column_stack
import seq_help
import gff_help

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


class SeqData(AbstractComputedData):

    AA_NUM_CODE = {aa:i for i, aa in enumerate(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 
                                       'G', 'H', 'I', 'L', 'K', 'M', 'F', 
                                       'P', 'S', 'T', 'W', 'Y', 'V'])}
    AA = sorted(AA_NUM_CODE, key=AA_NUM_CODE.get) # sort the key by value
    NUM_AA = 20

    def __init__(self, forceCompute=False): 
        rawInputName = 'MDB_Homo_sapiens/sequences.fasta'
        featureFileName = 'MDB_Homo_sapiens_seq_characters.tsv'
        super(SeqData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)


    def compute_scores(self):
        
        def score(index, pID):
            seqObj = seq_help.get_seq(index)
            seq = seqObj.seq
            startPos, endPos = 1, len(seq) # assume starting M removed
            chain = gff_help.get_longest_chain(pID) 
            if chain:
                startPos, endPos = chain

            nEndAA = seq[startPos]
            numLysine = seq.count('K')
            
            return (pID,
                    numLysine,
                #    nEndAA,
                ) + self.score_n_end(nEndAA)

        dtype=[('Uniprot', '|S20'),
               ('numLysine', int),
               #('nEndAA', '|S2'),
        ]
        
        dtype.extend(
            [(aa, int) for aa in self.AA]
        )

        allScores = [score(index, pid) for index,pid  in seq_help.get_all_prot()]  # the index is used to retrive the entry from the fasta file
        allScores = np.array(allScores,dtype=dtype)

        return allScores
    
    def score_n_end(self, nEndAA):
        c = self.AA_NUM_CODE.get(nEndAA,None)
        score = np.zeros(self.NUM_AA)
        if c != None:
            score[c] = 1
        return tuple(score)



class DisorderData(AbstractComputedData):

    def __init__(self, forceCompute=False): 
        rawInputName = 'MDB_Homo_sapiens/disorder_consensus.fasta'
        featureFileName = 'MDB_Homo_sapiens_disorder_consensus.csv'
        super(DisorderData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)


    def compute_scores(self):

        def term_disorder(aaDisorder, numAA, whichTerm):
            '''terminal disorder of this number of nterminal a.a.'''
            if whichTerm == 'C':
                terminal = aaDisorder[-numAA:]
            else:
                terminal = aaDisorder[:numAA]
            numDisorder = np.sum(terminal)
            return numDisorder / numAA

        def internal_disorder(aaDisorder, numAA):
            internalDisorderRegions = re.findall(r'1+',aaDisorder[40:])
            internaDisorderCnt = len([disOR for disOR in internalDisorderRegions if len(disOR)>numAA])
            return internaDisorderCnt

        varying_lengths = [30, 40, 50, 60, 70] 
        def score(record):
            aaDisorder = self.parse_disorder_consensus(record.seq) # e.g [1,1,1,0,...1,1], where 1 means disordered
            aaDisorderStr = ''.join(aaDisorder.astype(int).astype(str)) # e.g '1110...11'
            totalDisorderAA = np.sum(aaDisorder)
            totalDisorderRatio = totalDisorderAA / len(record.seq)

          #  import pdb; pdb.set_trace()
            
            ctd = [term_disorder(aaDisorder, l, 'C') for l in varying_lengths]
            ntd = [term_disorder(aaDisorder, l, 'N') for l in varying_lengths]
            intnld = [internal_disorder(aaDisorderStr, l) for l in varying_lengths]

            return (record.id, 
                    totalDisorderAA,
                    totalDisorderRatio,
                   ) + tuple(ntd) \
                   + tuple(ctd) \
                   + tuple(intnld)


        dtype=[('Uniprot', '|S20'),
               ('totalDisorderAA', int),
               ('totalDisorderRatio', float),
        ] + [('ntermDisorder'+str(l), float) for l in varying_lengths] \
        + [('ctermDisorder'+str(l), float) for l in varying_lengths] \
        + [('internaDisorderCnt'+str(l), float) for l in varying_lengths]


        with open(self.rawInputFilePath, "rU") as handle:
#            import pdb; pdb.set_trace()
            allScores = [score(rec) for rec in SeqIO.parse(handle, "fasta")]
            allScores = np.array(allScores,dtype=dtype)

        return allScores


    def parse_disorder_consensus(self, seq):
        # get array of amino acid -wise disorder consensus score
        aaDisorder = np.fromiter(seq, dtype=int, count=len(seq))
        return aaDisorder >= 5



class AbstractSitesData(AbstractComputedData):

    def setup(self):
        # not sure why genfromtxt refuses to use the first line as
        # header, when names=True, so doing this manually
        with open(self.rawInputFilePath,'r') as f:
            for x in range(3):
                f.readline()
            header = f.readline()
            header = header.rstrip().split('\t')
        self.sitesData = \
            np.genfromtxt(self.rawInputFilePath, 
                          dtype=None,
                          skip_header=3,
                          names=header,
                          delimiter='\t',
                        )

        self.column = \
            {cname:i for i,cname in enumerate(self.sitesData.dtype.names)}
  
        organism = self.sitesData['ORG']
        self.sitesData = self.sitesData[organism=='human']
        locate = {} # locate the entries by protein ID
        for i,entry in enumerate(self.sitesData):
            pid = self.getEntryID(entry)
            if not locate.get(pid):
                locate[pid] = []
            locate[pid].append(i)

        self.data_locate = locate
        return 

    def getEntryID(self, entry):
        return entry[self.column['ACC_ID']]

    def compute_scores(self):

        self.setup()

        def get__sites(pid):
            indexes = self.data_locate.get(pid)
            if indexes:
                sites = self.sitesData[indexes]
                sites['MOD_RSD']
            else:
                sites = []
            return sites

        def score(key, pID):
            sites = get__sites(pID)
            if sites == []:
                knownSites = 0
            else:
                knownSites = 1
            
            return (pID,
                len(sites),
                knownSites,
                )

        dtype=[('Uniprot', '|S20'),
               ('numSites', int),
               ('knownSites', int),
        ]
        #import pdb; pdb.set_trace() 
        # tag colunm names by the class
        dtype = [dtype[0]] + \
                [('%s_%s' %(self.class_abbr(),fn),t) for (fn,t) in dtype[1:]]
        

        pIDs = seq_help.get_all_prot()
        allScores = [score(key, pid) for key, pid in pIDs]
        allScores = np.array(allScores,dtype=dtype)

        return allScores

    def class_abbr(self):
        clssnm = self.__class__.__name__
        return clssnm.rstrip('SitesData')[:3]


class UbqSitesData(AbstractSitesData):

    def __init__(self, forceCompute=False): 
        rawInputName = 'Ubiquitination_site_dataset'
        featureFileName = 'Ubiquitination_site.tsv'
        self.abbr = 'Ubiq'
        super(UbqSitesData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)

class AcetSitesData(AbstractSitesData):

    def __init__(self, forceCompute=False): 
        rawInputName = 'Acetylation_site_dataset'
        featureFileName = 'Acetylation_site.tsv'
        super(AcetSitesData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)

class PhosphoSitesData(AbstractSitesData):

    def __init__(self, forceCompute=False): 
        rawInputName = 'Phosphorylation_site_dataset'
        featureFileName = 'Phosphorylation_site.tsv'
        super(PhosphoSitesData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)

class MetSitesData(AbstractSitesData):

    def __init__(self, forceCompute=False): 
        rawInputName = 'Methylation_site_dataset'
        featureFileName = 'Methylation_site.tsv'
        super(MetSitesData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)

        
if __name__ == '__main__':
    DisorderData(forceCompute=True)
    SeqData(forceCompute=True)
    UbqSitesData(forceCompute=True)
    MetSitesData(forceCompute=True)
    AcetSitesData(forceCompute=True)
    PhosphoSitesData(forceCompute=True)
                     
                 

