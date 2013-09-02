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



class Species(object):

    ALL_SPECIES = {}

    def __init__(self, speciesName, relPath=None):
        '''relPath: directory for species specific data, relative the RAW_INPUT_PATH'''
        speObj = self.ALL_SPECIES.get(speciesName)
        if speObj:
            raise Exception('Same Species already exists %s' % speObj)
        else:
            self.speciesName = speciesName
            self.input_path = relPath or speciesName
            self.dataClasses = []
            self.ALL_SPECIES[speciesName] = self
            
        super(Species, self).__init__()

    def register(self, dataClass):
        self.dataClasses.append(dataClass)

    def path(self, fileName):
        '''species specific data file path'''
        return os.path.join(self.input_path,
                            fileName)

    def __str__(self):
        return 'Species: %s @ %s' % (self.speciesName, self.path(''))
    


class SpeciesConnectorMetaclass(type): 
    def __new__(cls, clsname, bases, dct):
        featureClass = super(SpeciesConnectorMetaclass, cls).__new__(
            cls, clsname, bases, dct)
        if clsname != 'SPC':
            species = featureClass.SPECIES
            if not species or not isinstance(species, Species):
                raise Exception('Species Association Expected.') 
            species.register(featureClass)
        return featureClass



class SPC(object):
    __metaclass__ = SpeciesConnectorMetaclass
    SPECIES = None

    def path(self, fileName):
        return self.SPECIES.path(fileName)

HUMAN_SPECIES = Species('human')
class MyDisorderData(DisorderData, SPC):
    SPECIES = HUMAN_SPECIES

    def __init__(self, forceCompute=False):
        import pdb; pdb.set_trace()
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

def map_features(species):
    for dc in species.dataClasses:
        remap(dc())
    
    
if __name__ == '__main__':


    while(True):
        try:
            
            species_list = Species.ALL_SPECIES.values()
            t1 = lookup_table(species_list)
            print_table(t1, lambda item: item.speciesName)
            spCode = raw_input('Enter Species Code:\n')
            species = t1[int(spCode)]
        
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
