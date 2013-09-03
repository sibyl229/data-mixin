#from phl.settings import *
from settings import *
from compute_data import AbstractComputedData
from data import AbstractData


class Species(object):

    ALL_SPECIES = {}

    def __init__(self, moduleName, speciesName=None, relPath=None):
        '''
moduleName: where features for the species is defined. Usually under species/ subpackage
speciesName: a Human readable alias for your species
relPath: directory for species specific data, relative the RAW_INPUT_PATH'''
        speObj = self.ALL_SPECIES.get(moduleName)
        if speObj:
            raise Exception('Same Species already exists %s' % speObj)
        else:
            self.input_path = relPath or speciesName
            self.dataClasses = []
            self.ALL_SPECIES[moduleName] = self
            self.moduleName = moduleName
            self.speciesName = speciesName or moduleName.split('.')[-1]
        super(Species, self).__init__()

    def register(self, dataClass):
        self.dataClasses.append(dataClass)

    def path(self, fileName):
        '''species specific data file path'''
        return os.path.join(self.input_path,
                            fileName)

    def __str__(self):
        return 'Species: %s @ %s' % (self.speciesName, self.path(''))
    


class FeatureRegisterMetaclass(type): 
    def __new__(cls, clsname, bases, dct):
        featureClass = super(FeatureRegisterMetaclass, cls).__new__(
            cls, clsname, bases, dct)
        if issubclass(featureClass, AbstractComputedData) or \
           issubclass(featureClass, AbstractData):
            # tr
            species = featureClass.SPECIES
            if not species or not isinstance(species, Species):
                raise Exception('Species Association Expected.') 
            species.register(featureClass)
        return featureClass


class SpeciesConnectivityMixin(object):

    __metaclass__ = FeatureRegisterMetaclass
    
    SPECIES = None

    def path(self, fileName):
        return self.SPECIES.path(fileName)


