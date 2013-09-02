import sys
import os
from compute_data import *
from data import *
from abc import ABCMeta, abstractmethod
from DataRemap import remap

decode = {
    'q': sys.exit,
    #'c': compute_feature(sp)
}

SPECIES_BANK = []


class SpeciesBankMetaclass(type): 
    def __new__(cls, clsname, bases, dct):
        speciesClass = super(SpeciesBankMetaclass, cls).__new__(
            cls, clsname, bases, dct)
        if clsname != 'SpeciesFeaturesBase':
            if not speciesClass.speciesName:
                raise Exception('class variable speciesName is not defined in subclass')
            else:
                SPECIES_BANK.append(speciesClass)
                speciesClass.SPECIES_CLASS = cls
        return speciesClass



class SpeciesBase(object):

    __metaclass__ = SpeciesBankMetaclass

    def __init__(self):
        self.dataClasses = []
        self.setup()
        super(SpeciesFeaturesBase, self).__init__()

    def register(self, dataClass):
        self.dataClasses.append(dataClass)


    def setup(self):
        '''code for registering data classes'''
        pass


    def path(self, fileName):
        '''return path of the species folder relative to RAW_INPUT_PATH'''
        pass


class SpeciesConnectorMetaclass(type): 
    def __new__(cls, clsname, bases, dct):
        featureClass = super(SpeciesConnectorMetaclass, cls).__new__(
            cls, clsname, bases, dct)
        if clsname != 'SpeciesFeaturesBase':
            species = featureClass.SPECIES_CLASS
            if not species:
                raise Exception('Species Association Expected.') 
        return speciesClass


class HumanFeatures(SpeciesBase):
    speciesName = 'human'

    def path(self, fileName):
        '''return full path'''
        os.path.join('human', fileName)


class SPC(object):
    __metaclass__ = SpeciesConnectorMetaclass
        
class MyDisorderData(DisorderData, SPC):
    SPECIES_CLASS = HumanFeatures

    def __init__(self, forceCompute=False):
        super(MyDisorderData, self).__init__(
            rawInputName=self.path('MDB_Homo_sapiens/disorder_consensus.fasta'),
            featureFileName=self.path('disorder.tsv'),
            forceCompute=forceCompute)

    class MyDegrMotifData(DegrMotifData):    
        def __init__(self, forceCompute=False):
            super(MyDegrMotifData, self).__init__(
                rawInputName=self.path('MDB_Homo_sapiens/disorder_consensus.fasta'),
                featureFileName=self.path('degr_motif.tsv'),
                forceCompute=forceCompute)

    class MySeqData(SeqData):
        def __init__(self, forceCompute=False):
            super(MySeqData, self).__init__(
                rawInputName=self.path('MDB_Homo_sapiens/sequences.fasta'),
                featureFileName=self.path('seq_characters.tsv'),
                forceCompute=forceCompute)

    class MySecStructData(SecStructData):
        def __init__(self, forceCompute=False):
            super(MySecStructData, self).__init__(
                rawInputName=self.path('MDB_Homo_sapiens/annotations.fasta'),
                featureFileName=self.path('sec_struct.tsv'),
                forceCompute=forceCompute)

    class MyUbqSitesData(UbqSitesData):
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

    class MyPhosphoSitesData(PhosphoSitesData):
        def __init__(self, forceCompute=False):
            super(MyPhosphoSitesData, self).__init__(
                rawInputName = 'Phosphorylation_site_dataset',
                featureFileName = self.path('Phosphorylation_site.tsv'),
                organism='human',
                forceCompute=forceCompute
            )


    class MyMetSitesData(MetSitesData):
        def __init__(self, forceCompute=False):
            super(MyMetSitesData, self).__init__(
                rawInputName = 'Methylation_site_dataset',
                featureFileName = self.path('Methylation_site.tsv'),
                organism='human',
                forceCompute=forceCompute
            )

    def setup(self):
        featureDataClasses = [
            self.MyDegrMotifData,
            self.MySecStructData,
            self.MyDisorderData,
            self.MySeqData,
            self.MyUbqSitesData,
            self.MyMetSitesData,
            self.MyAcetSitesData,
            self.MyPhosphoSitesData,
        ]
        for c in featureDataClasses:
            self.register(c)
        return 
    
def lookup_table(lst):
        return {i:sp for i,sp in enumerate(lst)}

def print_table(dtable, to_string=None):
    for key, item in dtable.items():
        if to_string:
            itemStr = to_string(item)
        else:
            itemStr = str(item)
        print '%3d\t%s' % (key, itemStr)

def compute_features(species):

    t = lookup_table(species.dataClasses)
    print_table(t)
    response = raw_input(
'''Enter Data Code separated by white space:
Enter 'A' for all features and 'B' for Returning to Previous Step''')
    if re.search(r'b|B', response):
        return
    elif re.search(r'a|A',response):
        dataClasses = t.items()
    else:    
        dataIds = response.split()
        dataClasses = [t[int(i)] for i in dataIds]
    for dc in dataClasses:
        dc(forceCompute=True)
    return

def map_features(speciesClass):
    for dc in speciesClass.dataClasses:
        remap(dc())
    
    
if __name__ == '__main__':


    while(True):
        try:
            #import pdb; pdb.set_trace()
            t1 = lookup_table(SPECIES_BANK)
            print_table(t1, lambda item: item.speciesName)
            spCode = raw_input('Enter Species Code:\n')
            speciesClass = t1[int(spCode)]
            species = speciesClass()
        
            while(True):
                try:
                    opDict = {
                        'Compute Features': (lambda: compute_features(species)),
                        'Map Features': (lambda: map_features(species)),
                        #   'Combine Features': lambda: ,
                    }
                    t2 = lookup_table(opDict.keys())
                    #import pdb; pdb.set_trace()

                    print_table(t2)
                    opCode = raw_input('Enter Operation Code:\n')
                    opName = t2[int(opCode)]
                    opDict[opName]()
    
                except KeyError, ValueError:
                    print 'Invalid Operation'
                
        except KeyboardInterrupt:
            sys.exit()
        except KeyError, ValueError:
            print 'Invalid Species Code!'
