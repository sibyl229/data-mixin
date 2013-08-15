from __future__ import division
import os
import numpy as np
import re
from Bio import SeqIO
from abc import ABCMeta, abstractmethod
from settings import *
from helpers import column_stack

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
        self.proteinGffPath = os.path.join(RAW_INPUT_PATH, 'human_uniprot.gff.tsv')
        self.gffHeader = ['Uniprot','source','feature','start','end','attr']
        self.setup_gf()
        super(SeqData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)


    def compute_scores(self):
        
        def is_refseq(fasta_header): 
            match = re.search(r'(\w+)\|refseq', fasta_header)
            if match:
                refseqId = match.group(1)
            else:
                refseqId = 'non-refseq'
            return refseqId
        
        record_dict = SeqIO.index(self.rawInputFilePath, "fasta")
        
        def score(refseqIndex):
            pID = is_refseq(refseqIndex)
            seq = record_dict[refseqIndex].seq
            
            startPos, endPos = 0, len(seq)
            located = self.get_longest_chain(pID)
            if located:
                startPos, endPos = located 

            nEndAA = seq[startPos]
            numLysine = seq.count('K')
            
            return (pID,
                    numLysine,
                #    nEndAA,
                ) + self.score_n_end(nEndAA)

        dtype=[('Uniprot', '|S20'),
               ('numLysine', int),
               #('nEndAA', '|S2'),
               # ('totalDisorderRatio', float),
               # ('ntermDisorder', float),
               # ('internaDisorderCnt', int),
        ]
        
        dtype.extend(
            [(aa, int) for aa in self.AA]
        )

        get_seq_id = np.vectorize(is_refseq)
        seqIndexes = np.array(record_dict.keys())
        refseqIndexes = seqIndexes[get_seq_id(seqIndexes) != 'non-refseq']
        allScores = [score(rIndex) for rIndex in refseqIndexes]
        allScores = np.array(allScores,dtype=dtype)
#        import pdb; pdb.set_trace() 

        return allScores

    
    def get_longest_chain(self, protId):
        protGf = self.get_protein_gf(protId)
        isChain = protGf['feature']=='Chain'
        protGf = protGf[isChain]
        if len(protGf) > 0:
            chainLengths = protGf['end'] - protGf['start']
            indexLongest = np.argmax(chainLengths)
            start, end = protGf[['start', 'end']][indexLongest]

            # change 1-indexed position to 0-indexed position
            # AND change closed range [start, end] to open range [start, end)
            # so they follow programming conventions
            start, end = start-1, end
            return start, end
        else:
            return None


    def setup_gf(self):
        gff = np.genfromtxt(self.proteinGffPath, dtype=None, 
                      usecols=[0,1,2,3,4,8],
                      missing='.', delimiter='\t',
                      names=self.gffHeader)
        locate = {}
        for i,entry in enumerate(gff):
            pid = entry[0]
            if not locate.get(pid):
                locate[pid] = []
            locate[pid].append(i)

        self.gff = gff
        self.gff_locate = locate
        return gff, locate

    
    def get_protein_gf(self, protId):
        indices = self.gff_locate.get(protId, [])
        return self.gff[indices]
            

    
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


class UbqSitesData(AbstractComputedData):


    def __init__(self, forceCompute=False): 
        rawInputName = 'Ubiquitination_site_dataset'
        featureFileName = 'Ubiquitination_site.tsv'
        super(UbqSitesData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)


    def setup(self):
        # not sure why genfromtxt refuses to use the first line as
        # header, when names=True, so doing this manually
        with open(self.rawInputFilePath,'r') as f:
            header = f.readline()
            header = header.rstrip().split('\t')
        self.ubqData = \
            np.genfromtxt(self.rawInputFilePath, 
                          dtype=None, 
                          names=header,
                          delimiter='\t',
                        )

        self.column = \
            {cname:i for i,cname in enumerate(self.ubqData.dtype.names)}
  
        organism = self.ubqData['ORG']
        self.ubqData = self.ubqData[organism=='human']
        locate = {} # locate the entries by protein ID
        for i,entry in enumerate(self.ubqData):
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

        def get_ubq_sites(pid):
            indexes = self.data_locate.get(pid)
            if indexes:
                sites = self.ubqData[indexes]
                sites['MOD_RSD']
            else:
                sites = []
            return sites

        def score(pID):
            ubqsites = get_ubq_sites(pID)
            if ubqsites == []:
                knownUbqSites = 0
            else:
                knownUbqSites = 1
            
            return (pID,
                len(ubqsites),
                knownUbqSites,
                )

        dtype=[('Uniprot', '|S20'),
               ('numUbqSites', int),
               ('knownUbqSites', int),
               # ('totalDisorderRatio', float),
               # ('ntermDisorder', float),
               # ('internaDisorderCnt', int),
        ]
        
#        import pdb; pdb.set_trace() 
        pIDs = set([self.getEntryID(entry) for entry in self.ubqData])
        allScores = [score(pid) for pid in pIDs]
        allScores = np.array(allScores,dtype=dtype)

        return allScores

    
    def get_protein_(self, protId):
        indices = self.gff_locate.get(protId, [])
        return self.ubqData[indices]
    

        
if __name__ == '__main__':
#    SeqData(forceCompute=True)
    UbqSitesData(forceCompute=True)
    DisorderData(forceCompute=True)

