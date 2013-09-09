import os
from Bio import SeqIO
from compute_data import *
from data import *
from species.base import *
from helpers import column_stack


class UniprotToSGD(object):
    '''maps Uniprot/Swissprot id to SGD id'''
    def __init__(self, species, fromIDType='UniProt/Swiss-Prot ID'):
        self.SPECIES = species
        idMappingFilePath = os.path.join(RAW_INPUT_PATH,
                                         self.SPECIES.rel_path('dbxref.tab'))
        header = [
            'DBXREF_ID',
            'DBXREF_ID_source',
            'DBXREF_ID_type',
            'Systematic_Name',
            'SGDID',
            'Gene_Name',
        ]

        mapping = np.genfromtxt(idMappingFilePath, 
                                delimiter='\t', 
                                dtype=None,
                                names=header)
        relevantRows = (mapping['DBXREF_ID_type']==fromIDType)
        cleanMapping = mapping[relevantRows][['DBXREF_ID', 'Systematic_Name']]
        self.idLookup = {uniprotID: sysName for uniprotID, sysName in cleanMapping}
         
    
    def translateID(self, anID):
        return self.idLookup.get(anID, '')
        
    def mapIDs(self, ids):
        mappedIDs = [self.translateID(idx) for idx in ids]
        return np.array(mappedIDs)

    def annotate(self, data, currIDcolName):
        currIDs = data[currIDcolName]
        idColumns = np.zeros((data.shape[0],), 
                             dtype=[('PrimaryID','|S9'), ('OtherID','|S28')])
        mappedIDs = self.mapIDs(currIDs)
        idColumns['PrimaryID'] = mappedIDs
        idColumns['OtherID'] = currIDs
        otherColNames = [nm for nm in data.dtype.names if nm != currIDcolName]
        annotated = column_stack(idColumns,
                                 data[otherColNames])
        return annotated[ mappedIDs!='' ]


class MySPC(SpeciesConnectivityMixin):
    SPECIES = Species.ALL_SPECIES[__name__]
    ### Register a few context for the species
    SPECIES.context.update({
        'gffName': '',
        'goaName': 'gene_association', # gene ontology
        'taxonNum': '4932',
        'seqName': 'orf_trans.fasta', 
    })

## class that contains Half Lifes
## used as template for organizing other features
class DegRateData(AbstractData, MySPC):
    CONTAINS_HALF_LIFE_DATA = True
    HALF_LIFE_COL_NAME = 'Corrected_Half_Life'
    def __init__(self):
        originalFileName = self.path('pnas_0605420103_SuppDataSet.txt')
        featureFileName = self.path('pnas_0605420103_SuppDataSet_cleaned.txt')
        self.clean_file(originalFileName, featureFileName) 
        idColName = 'ORF'
        super(DegRateData, self).__init__(
            featureFileName, idColName, usecols=[0,3])

    def clean_file(self, originalFileName, cleanedFileName):
        '''fix some weirdness of the file, such as inconsistent column number in some entries'''
        originalFilePath = os.path.join(RAW_INPUT_PATH, originalFileName)
        cleanedFilePath = os.path.join(RAW_INPUT_PATH, cleanedFileName)
        with open(originalFilePath, 'r') as infile:
            lines = infile.readlines()
            with open(cleanedFilePath, 'w') as outfile:
                for i, ln in enumerate(lines[5:]):
                    outfile.write(ln.rstrip()+'\r\n')
        


#Below are other feature classes
class MyDisorderData(DisorderData, MySPC):

    def __init__(self, forceCompute=False):
        originalRawInputName = self.path('disorder_consensus.fasta')
        rawInputName = self.path('disorder_consensus_cleaned.fasta')
        self.clean_file(originalRawInputName, rawInputName)
        super(MyDisorderData, self).__init__(
            rawInputName=rawInputName,
            featureFileName=self.path('disorder.tsv'),
            forceCompute=forceCompute)

    def clean_file(self, originalFileName, cleanedFileName):
        '''fix some weirdness of the file, such as inconsistent column number in some entries'''
        originalFilePath = os.path.join(RAW_INPUT_PATH, originalFileName)
        cleanedFilePath = os.path.join(RAW_INPUT_PATH, cleanedFileName)
        idMapper = UniprotToSGD(self.SPECIES)
        uniqueId = set()
        with open(cleanedFilePath, 'w') as outfile:
            for record in SeqIO.parse(open(originalFilePath, 'r'), 'fasta') :
                uniprotID = record.id
                sysName = idMapper.translateID(uniprotID)
                if sysName and not sysName in uniqueId:
                    uniqueId.add(sysName)
                    seq = str(record.seq)
                    newRecoreStr = '>%s\n%s\n' % (sysName, seq)
                    outfile.write(newRecoreStr)
        return

    # def compute_scores(self, scores, backupPath):
    #     scoresNew = UniprotToSGD(self.SPECIES).annotate(scores, 'Uniprot')
    #     return scoresNew


class MyDegrMotifData(DegrMotifData, MySPC):    
    def __init__(self, forceCompute=False):
        super(MyDegrMotifData, self).__init__(
            rawInputName=self.path('disorder_consensus_cleaned.fasta'),
            featureFileName=self.path('degr_motif.tsv'),
            forceCompute=forceCompute)

# class MySeqData(SeqData, MySPC):
#     def __init__(self, forceCompute=False):
#         super(MySeqData, self).__init__(
#             rawInputName=self.path('MDB_Homo_sapiens/sequences.fasta'),
#             featureFileName=self.path('seq_characters.tsv'),
#             forceCompute=forceCompute)

# class MySecStructData(SecStructData, MySPC):
#     def __init__(self, forceCompute=False):
#         super(MySecStructData, self).__init__(
#             rawInputName=self.path('MDB_Saccharomyces_cerevisiae_(strain_ATCC_204508_|_S288c)/annotations.fasta'),
#             featureFileName=self.path('sec_struct.tsv'),
#             forceCompute=forceCompute)

class PropertiesData(AbstractData, MySPC):
    headerInfo = \
'''
0 	Orf_name
1 	SGDID
2 	MOLECULAR WEIGHT (in Daltons)
3 	PI
4 	CAI (Codon Adaptation Index)
5 	PROTEIN LENGTH
6 	N TERM SEQ
7 	C TERM SEQ
8 	CODON BIAS
9 	ALA
10 	ARG
11 	ASN
12 	ASP
13 	CYS
14 	GLN
15 	GLU
16 	GLY
17 	HIS
18 	ILE
19 	LEU
20 	LYS
21 	MET
22 	PHE
23 	PRO
24 	SER
25 	THR
26 	TRP
27 	TYR
28 	VAL
29 	FOP SCORE (Frequency of Optimal Codons)
30 	GRAVY SCORE (Hydropathicity of Protein)
31 	AROMATICITY SCORE (Frequency of aromatic amino acids: Phe, Tyr, Trp)
32 	Feature type (ORF classification: Verified, Uncharacterized, Dubious'''

    headerAll = headerInfo.strip().split('\n')

    def __init__(self, forceCompute=False):
        featureFileName = self.path('protein_properties.tab')
        idColName = 'Orf_name'
        usecols = sorted(list(set(range(33)) - set([1, 6, 7, 32])))
        useColNames = [self.abbr(self.headerAll[col]) \
                       for col in usecols]
        super(PropertiesData, self).__init__(
            featureFileName, idColName, usecols=usecols, names=useColNames)

    def abbr(self, name):
        '''make shorter field name'''
        print name
        name = re.search(r'\d+\s+(.*)', name).group(1)
        name = re.sub(r'\(.*\)', '', name)
        name = re.sub('\s', '_', name.strip())
        return name
    
    @staticmethod
    def is_full_entry(t):
        '''check if a tuple has nan value'''
        t = tuple(t)
        return not np.isnan(np.sum(t))

    def _extract_data(self): 
        '''remove records with missing values'''
        idCol, data = super(PropertiesData, self)._extract_data()
        fullEntries = np.array([self.is_full_entry(row) for row in data])

        self.rawdata = self.rawdata[fullEntries]
        return super(PropertiesData, self)._extract_data()


