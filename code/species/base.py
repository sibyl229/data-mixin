#from phl.settings import *
from settings import *
from compute_data import AbstractComputedData
from data import AbstractData
import os


class Species(object):

    ALL_SPECIES = {}

    def __init__(self, moduleName, speciesName=None, relPath=None):
        '''
moduleName: where features for the species is defined. Usually under species/ subpackage
speciesName: a Human readable alias for your species
relPath: directory for species specific data, relative to the RAW_INPUT_PATH'''
        speObj = self.ALL_SPECIES.get(moduleName)
        if speObj:
            raise Exception('Same Species already exists %s' % speObj)
        else:
            self.relPath = relPath or speciesName
            self.dataClasses = []
            self.HLDataClass = None
            self.ALL_SPECIES[moduleName] = self
            self.moduleName = moduleName
            self.speciesName = speciesName or moduleName.split('.')[-1]
            # create directory to store species specific cleaned feature
            # files
            make_sure_path_exists(os.path.join(CLEAN_INPUT_PATH, 
                                               self.rel_path()))
        super(Species, self).__init__()

    def register(self, dataClass):
        self.dataClasses.append(dataClass)

    def get_computed_data_classes(self):
        return [dc for dc in self.dataClasses \
            if issubclass(dc, AbstractComputedData)]

    def rel_path(self, fileName=''):
        '''species specific data file path'''
        return os.path.join(self.relPath,
                            fileName)


    def __str__(self):
        return 'Species: %s @ %s' % (self.speciesName, self.path(''))
    


class FeatureRegisterMetaclass(type): 
    def __new__(cls, clsname, bases, dct):
        featureClass = super(FeatureRegisterMetaclass, cls).__new__(
            cls, clsname, bases, dct)
        if issubclass(featureClass, AbstractComputedData) or \
           issubclass(featureClass, AbstractData):
            #
            species = featureClass.SPECIES
            if not species or not isinstance(species, Species):
                raise Exception('Species Association Expected.') 
            species.register(featureClass)

            if getattr(featureClass, 
                       'CONTAINS_HALF_LIFE_DATA',
                       None):
                if species.HLDataClass:
                    raise Exception('More than one feature data object tagged with CONTAINS_HALF_LIFE_DATA flag.')
                species.HLDataClass = featureClass
        return featureClass


class SpeciesConnectivityMixin(object):

    __metaclass__ = FeatureRegisterMetaclass
    
    SPECIES = None

    def path(self, fileName):
        return self.SPECIES.rel_path(fileName)


