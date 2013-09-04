import os
from compute_data import *
from data import *
from species.base import *



class MySPC(SpeciesConnectivityMixin):
    SPECIES = Species.ALL_SPECIES[__name__]
    ### Register a few context for the species
    SPECIES.context.update({
        'gffName': 'uniprot-taxonomy%3A9606+AND+reviewed%3Ayes+AND+keyword%3A181.gff',
        'fastaName': '',
    })

## class that contains Half Lifes
## used as template for organizing other features
class DegRateData(AbstractData, MySPC):
    CONTAINS_HALF_LIFE_DATA = True
    HALF_LIFE_COL_NAME = 'halflife_t12_in_h'
    def __init__(self):
        featureFileName = self.path('pr101183k_si_002_HeLa.csv')
        idColName = 'Uniprot'
        usecols = [3,11]
        super(DegRateData, self).__init__(
            featureFileName, idColName, usecols=usecols)


## Below are other feature classes
class MyDisorderData(DisorderData, MySPC):

    def __init__(self, forceCompute=False):
        super(MyDisorderData, self).__init__(
            rawInputName=self.path('MDB_Homo_sapiens/disorder_consensus.fasta'),
            featureFileName=self.path('disorder.tsv'),
            forceCompute=forceCompute)

class MyDegrMotifData(DegrMotifData, MySPC):    
    def __init__(self, forceCompute=False):
        super(MyDegrMotifData, self).__init__(
            rawInputName=self.path('MDB_Homo_sapiens/disorder_consensus.fasta'),
            featureFileName=self.path('degr_motif.tsv'),
            forceCompute=forceCompute)

class MySeqData(SeqData, MySPC):
    def __init__(self, forceCompute=False):
        super(MySeqData, self).__init__(
            rawInputName=self.path('MDB_Homo_sapiens/sequences.fasta'),
            featureFileName=self.path('seq_characters.tsv'),
            forceCompute=forceCompute)

class MySecStructData(SecStructData, MySPC):
    def __init__(self, forceCompute=False):
        super(MySecStructData, self).__init__(
            rawInputName=self.path('MDB_Homo_sapiens/annotations.fasta'),
            featureFileName=self.path('sec_struct.tsv'),
            forceCompute=forceCompute)

class MyUbqSitesData(UbqSitesData, MySPC):
    def __init__(self, forceCompute=False):
        super(MyUbqSitesData, self).__init__(            
            rawInputName = 'Ubiquitination_site_dataset', #input shared among species
            featureFileName = self.path('Ubiquitination_site.tsv'),
            organism='human',
            forceCompute=forceCompute)

class MyAcetSitesData(AcetSitesData):
    def __init__(self, forceCompute=False):
        super(MyAcetSitesData, self).__init__(
            rawInputName = 'Acetylation_site_dataset',
            featureFileName = self.path('Acetylation_site.tsv'), 
            organism='human',
            forceCompute=forceCompute
        )

class MyPhosphoSitesData(PhosphoSitesData, MySPC):
    def __init__(self, forceCompute=False):
        super(MyPhosphoSitesData, self).__init__(
            rawInputName = 'Phosphorylation_site_dataset',
            featureFileName = self.path('Phosphorylation_site.tsv'),
            organism='human',
            forceCompute=forceCompute
        )

class MyMetSitesData(MetSitesData, MySPC):
    def __init__(self, forceCompute=False):
        super(MyMetSitesData, self).__init__(
            rawInputName = 'Methylation_site_dataset',
            featureFileName = self.path('Methylation_site.tsv'),
            organism='human',
            forceCompute=forceCompute
        )


