from __future__ import division
import os
import numpy as np
import re
from Bio import SeqIO
from abc import ABCMeta, abstractmethod
from settings import *
from gff_help import Gff
from seq_help import MySequence

class AbstractComputedData(object):
    '''Similar to AbstractData class, AbstractComputedData also implements get_all_data() used by remap().
Different from AbstractData class, this one is extended by classes that compute feature scores 
(and backs them up in a feature file).
 '''    
    
    @abstractmethod
    def __init__(self, inputFileName, featureFileName, forceCompute=False):
        self.rawInputFilePath = os.path.join(RAW_INPUT_PATH, inputFileName)
        self.featureFileName = featureFileName
        self.featureFilePath = os.path.abspath(
            os.path.join(RAW_INPUT_PATH, 
                         self.featureFileName))
 
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
        print 'computing feature: ' + self.__class__.__name__ + '...'
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

    def __init__(self, forceCompute=False,
        rawInputName = 'MDB_Homo_sapiens/sequences.fasta',
        featureFileName = 'MDB_Homo_sapiens_seq_characters.tsv'):
        super(SeqData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)


    def compute_scores(self):

        gff = Gff(self.SPECIES) # general feature file
        seqHandle = MySequence(self.SPECIES)

        def score(index, pID):
            seqObj = seqHandle.get_seq(index)
            seq = seqObj.seq
            #startPos, endPos = 1, len(seq) # assume starting M removed
            startPos, endPos = 0, len(seq) # assume starting with M
            chain = gff.get_longest_chain(pID) 
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

        allScores = [score(index, pid) for index,pid  in seqHandle.get_all_prot()]  # the index is used to retrive the entry from the fasta file
        allScores = np.array(allScores,dtype=dtype)

        return allScores
    
    def score_n_end(self, nEndAA):
        c = self.AA_NUM_CODE.get(nEndAA,None)
        score = np.zeros(self.NUM_AA)
        if c != None:
            score[c] = 1
        return tuple(score)



class DisorderData(AbstractComputedData):

    def __init__(self, forceCompute=False, 
                 rawInputName='MDB_Homo_sapiens/disorder_consensus.fasta',
                 featureFileName='MDB_Homo_sapiens_disorder_consensus.csv'):
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
            aaDisorderStr = self.to_disorder_str(aaDisorder)
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

    @staticmethod
    def parse_disorder_consensus(seq):
        # get array of amino acid -wise disorder consensus score
        aaDisorder = np.fromiter(seq, dtype=int, count=len(seq))
        return aaDisorder >= 5

    @staticmethod
    def to_disorder_str(disorderArray):
        return \
            ''.join(disorderArray.astype(int).astype(str)) # e.g '1110...11'

class SecStructData(AbstractComputedData):

    def __init__(self, forceCompute=False,
                 rawInputName='MDB_Homo_sapiens/annotations.fasta',
                 featureFileName='MDB_Homo_sapiens_sec_struct.csv'):
        super(SecStructData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)


    def compute_scores(self):

        def score(record):
            seq = record.seq
            fracHelix = self.frac(seq, r'[HGI]')
            fracSheet = self.frac(seq, r'[BE]')
            fracBend = self.frac(seq, r'[STL]')
            fracLoop = self.frac(seq, r'C')
            return (fracHelix, fracSheet, fracBend, fracLoop)


        dtype=[('Uniprot', '|S20'),
               ('fracHelix', float),
               ('fracSheet', float),
               ('fracBend', float),
               ('fracLoop', float),
        ]

        fastaPath = self.rawInputFilePath
        recordDict = SeqIO.index(fastaPath, "fasta")
        proteinIds = set([])
        allScores = []
        for key in recordDict.keys():
            match = re.search(r'(\w+)\|.*dssp', key)
            if match:
                pid = match.group(1)
                if pid not in proteinIds:
                    seqRec = recordDict.get(key)
                    allScores.append((pid,) + score(seqRec))
                    proteinIds.add(pid)
        
        allScores = np.array(allScores, dtype=dtype)
        return allScores        
                
    @staticmethod
    def frac(seq, regex):
        marker = 'A'
        seqNew = re.sub(regex, marker, str(seq))
        count = seqNew.count(marker)
        return count / len(seq)


class DegrMotifData(DisorderData):

    def __init__(self, forceCompute=False,
                 rawInputName='MDB_Homo_sapiens/disorder_consensus.fasta',
                 featureFileName='MDB_Homo_sapiens_degr_motif.csv'): 
        AbstractComputedData.__init__(self,
            rawInputName, featureFileName, forceCompute=forceCompute)


    def compute_scores(self):
        seqHandle = MySequence(self.SPECIES)

        def score(disorderSeq, seq):
            dregions = self.disorder_regions(disorderSeq, seq)
            KENPattern = r'KEN...N'
            DPattern = r'R..L....N'
            kCount = 0
            dCount = 0
            for rgn in dregions:
                kCount += len(re.findall(KENPattern, rgn))
                dCount += len(re.findall(DPattern, rgn))
            
            # if kCount > 0 or dCount > 0:
            #     import pdb; pdb.set_trace()
            kCountWhole = len(re.findall(KENPattern, seq))
            dCountWhole = len(re.findall(DPattern, seq))
            
            return (kCount, dCount, kCountWhole, dCountWhole)


        dtype=[('Uniprot', '|S20'),
               ('numKENbox', int),
               ('numDbox', int),
               ('numKENboxWhole', int),
               ('numDboxWhole', int),
        ]

        allScores = []
        disorderFastaDict = SeqIO.index(self.rawInputFilePath, 'fasta')
        for key, pid in seqHandle.get_all_prot():
            disorderConsensus = disorderFastaDict.get(pid)
            if disorderConsensus:
                seqRec = seqHandle.get_seq(key)
                seq = str(seqRec.seq)
                allScores.append(
                    (pid,) + score(disorderConsensus, seq))

        allScores = np.array(allScores,dtype=dtype)
        return allScores

    @classmethod
    def disorder_regions(cls, disorderSeq, seq):
        da = cls.parse_disorder_consensus(disorderSeq)
        ds = cls.to_disorder_str(da)
        matches = re.finditer('1+', ds)
        regions = [seq[m.start(): m.end()] for m in matches]
        #import pdb; pdb.set_trace()
        return regions
        


class AbstractSitesData(AbstractComputedData):
    '''Deals with multispecies Ubiquitination, acetylation site data file'''
    
    @abstractmethod
    def __init__(self,
                 rawInputName,
                 featureFileName,
                 organism,
                 forceCompute=False): 
        featureFileName = re.sub(r'(.+\.)(\w+)$', 
                                 r'\1%s.\2' % organism, 
                                 featureFileName) # insert organism name before the file extenstion
        self.organism = organism # used to filter entries in multiorganism data
        super(AbstractSitesData, self).__init__(
             rawInputName, featureFileName, forceCompute=forceCompute)


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
        self.sitesData = self.sitesData[organism==self.organism]
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
        

        pIDs = seqHandle.get_all_prot()
        allScores = [score(key, pid) for key, pid in pIDs]
        allScores = np.array(allScores,dtype=dtype)

        return allScores

    def class_abbr(self):
        clssnm = self.__class__.__name__
        return clssnm.rstrip('SitesData')[:3]


class UbqSitesData(AbstractSitesData):

    def __init__(self, forceCompute=False,
                 rawInputName = 'Ubiquitination_site_dataset',
                 featureFileName = 'Ubiquitination_site.tsv',
                 organism='human'): 
        
        self.abbr = 'Ubiq'
        super(UbqSitesData, self).__init__(
             rawInputName, featureFileName, organism, forceCompute=forceCompute)

class AcetSitesData(AbstractSitesData):

    def __init__(self, forceCompute=False,
                 rawInputName = 'Acetylation_site_dataset',
                 featureFileName = 'Acetylation_site.tsv', 
                 organism='human'): 

        super(AcetSitesData, self).__init__(
             rawInputName, featureFileName, organism, forceCompute=forceCompute)

class PhosphoSitesData(AbstractSitesData):

    def __init__(self, forceCompute=False,
                 rawInputName = 'Phosphorylation_site_dataset',
                 featureFileName = 'Phosphorylation_site.tsv',
                 organism='human'): 

        super(PhosphoSitesData, self).__init__(
             rawInputName, featureFileName, organism, forceCompute=forceCompute)

class MetSitesData(AbstractSitesData):

    def __init__(self, forceCompute=False,
                 rawInputName = 'Methylation_site_dataset',
                 featureFileName = 'Methylation_site.tsv',
                 organism='human'): 

        super(MetSitesData, self).__init__(
             rawInputName, featureFileName, organism, forceCompute=forceCompute)
                     
                 


