import os
from compute_data import *
from data import *
from species.base import *
from helpers import column_stack



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
    HALF_LIFE_COL_NAME = 'Corrected Half Life'
    def __init__(self):
        originalFileName = self.path('pnas_0605420103_SuppDataSet.txt')
        featureFileName = self.path('pnas_0605420103_SuppDataSet_cleaned.txt')
        self.clean_file(originalFileName, featureFileName) 
        idColName = 'ORF'
        super(DegRateData, self).__init__(
            featureFileName, idColName)

    def clean_file(self, originalFileName, cleanedFileName):
        '''fix some weirdness of the file, such as inconsistent column number in some entries'''
        originalFilePath = os.path.join(RAW_INPUT_PATH, originalFileName)
        cleanedFilePath = os.path.join(RAW_INPUT_PATH, cleanedFileName)
        with open(originalFilePath, 'r') as infile:
            lines = infile.readlines()
            with open(cleanedFilePath, 'w') as outfile:
                for i, ln in enumerate(lines[5:]):
                    # lnlen = len(ln.split('\t'))
                    # if lnlen != 4:
                    #     print i, lnlen, ln.split('\t')
                    outfile.write(ln.rstrip()+'\r\n')
        

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
         
    
    def mapIDs(self, ids):
        mappedIDs = [self.idLookup.get(idx, '') for idx in ids]
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
        import pdb; pdb.set_trace()        
        return annotated[ mappedIDs!='' ]
                

        
                                    
        
        

#Below are other feature classes
class MyDisorderData(DisorderData, MySPC):

    def __init__(self, forceCompute=False):
        super(MyDisorderData, self).__init__(
            rawInputName=self.path('disorder_consensus.fasta'),
            featureFileName=self.path('disorder.tsv'),
            forceCompute=forceCompute)

    def backup_scores(self, scores, backupPath):
        scoresNew = UniprotToSGD(self.SPECIES).annotate(scores, 'Uniprot')
        import pdb; pdb.set_trace()        
        super(MyDisorderData, self).backup_scores(scoresNew, backupPath)

# class MyDegrMotifData(DegrMotifData, MySPC):    
#     def __init__(self, forceCompute=False):
#         super(MyDegrMotifData, self).__init__(
#             rawInputName=self.path('MDB_Homo_sapiens/disorder_consensus.fasta'),
#             featureFileName=self.path('degr_motif.tsv'),
#             forceCompute=forceCompute)

# class MySeqData(SeqData, MySPC):
#     def __init__(self, forceCompute=False):
#         super(MySeqData, self).__init__(
#             rawInputName=self.path('MDB_Homo_sapiens/sequences.fasta'),
#             featureFileName=self.path('seq_characters.tsv'),
#             forceCompute=forceCompute)

# class MySecStructData(SecStructData, MySPC):
#     def __init__(self, forceCompute=False):
#         super(MySecStructData, self).__init__(
#             rawInputName=self.path('MDB_Homo_sapiens/annotations.fasta'),
#             featureFileName=self.path('sec_struct.tsv'),
#             forceCompute=forceCompute)

# class MyUbqSitesData(UbqSitesData, MySPC):
#     def __init__(self, forceCompute=False):
#         super(MyUbqSitesData, self).__init__(            
#             rawInputName = 'Ubiquitination_site_dataset', #input shared among species
#             featureFileName = self.path('Ubiquitination_site.tsv'),
#             organism='human',
#             forceCompute=forceCompute)

# class MyAcetSitesData(AcetSitesData):
#     def __init__(self, forceCompute=False):
#         super(MyAcetSitesData, self).__init__(
#             rawInputName = 'Acetylation_site_dataset',
#             featureFileName = self.path('Acetylation_site.tsv'), 
#             organism='human',
#             forceCompute=forceCompute
#         )

# class MyPhosphoSitesData(PhosphoSitesData, MySPC):
#     def __init__(self, forceCompute=False):
#         super(MyPhosphoSitesData, self).__init__(
#             rawInputName = 'Phosphorylation_site_dataset',
#             featureFileName = self.path('Phosphorylation_site.tsv'),
#             organism='human',
#             forceCompute=forceCompute
#         )

# class MyMetSitesData(MetSitesData, MySPC):
#     def __init__(self, forceCompute=False):
#         super(MyMetSitesData, self).__init__(
#             rawInputName = 'Methylation_site_dataset',
#             featureFileName = self.path('Methylation_site.tsv'),
#             organism='human',
#             forceCompute=forceCompute
#         )


