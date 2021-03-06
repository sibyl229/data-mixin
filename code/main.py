import sys
import os
from compute_data import *
from data import *
#from abc import ABCMeta, abstractmethod
from DataRemap import remap
from species.base import Species
from learning_setup import combine_features
from linear import simple_linear
from settings import *


####
# The Command Line User Interface
# 
# including the selection of a species and:
#  0	Compute Features
#  1	Map Features
#  2	Combine Features
#  3	Predict
####

### Register species modules HERE!
### In the format of Species(modulename, alias), 
###where alias is also the species specific subdirectory in raw_inputs/
### OR Species(modulename, alias, directory), to specifiy another directory
###base on its relative path from raw_inputs/
### A module specifies data/features for a species, including: 
###their file locations and preprocessing steps
species_choices = [
   # Species('species.human', 'Human'), # unless you get data from me, too big for github
    Species('species.yeast', 'Yeast'),
]


### helpers for command line interface    
def lookup_table(lst, f=lambda item: item):
    return {i:f(sp) for i,sp in enumerate(lst)}

def print_table(dtable, to_string=None):
    for key, item in dtable.items():
        if to_string:
            itemStr = to_string(item)
        else:
            itemStr = str(item)
        print '%3d\t%s' % (key, itemStr)
    print 'Ctrl-c\tExit'

def parse_int(inStr):
    matchInt = re.search(r'(\d+)',inStr)
    intVal = None
    if matchInt:
        intVal = int(matchInt.group(1))
    return intVal
        

def compute_features(species):

    t = lookup_table(species.get_computed_data_classes())
    print_table(t)
    response = raw_input(
        '''Enter Data Code separated by white space:
Enter 'A' for all features and 'B' for Returning to Previous Step
''')
    if re.search(r'b|B', response):
        return
    elif re.search(r'a|A',response):
        dataClasses = t.values()
    else:    
        dataIds = response.split()
        dataIds = [parse_int(i) for i in dataIds if parse_int(i) != None]
        dataClasses = [t[i] for i in dataIds]
    for dc in dataClasses:
        dc(forceCompute=True)
    return

def predict(species):
    normalizedFeaturePath = \
        os.path.join(CLEAN_INPUT_PATH,
                     species.rel_path('normalized_features.tsv'))
    simple_linear(normalizedFeaturePath)
    return


class MenuReturn(Exception):
    pass

def raise_menu_return_exception():
    raise MenuReturn()

### command line interface
    
if __name__ == '__main__':


    while(True):
        try:
            
            
            t1 = lookup_table(species_choices)
            print_table(t1, lambda item: item.speciesName)
            spCode = raw_input('Enter Species Code:\n')
            species = t1[parse_int(spCode)]
            species_module = __import__(species.moduleName)

        
            while(True):
                try:
                    operations = [
                        ('Compute Features', lambda: compute_features(species)),
                        ('Map Features', lambda: remap(species)),
                        ('Combine Features', lambda: combine_features(species)),
                        ('Predict', lambda: predict(species)),
                        ('Back to Previous', raise_menu_return_exception)
                    ]
                    t2 = lookup_table(operations, f=lambda itm: itm[0])
                    #import pdb; pdb.set_trace()

                    print_table(t2)
                    opCode = raw_input('Enter Operation Code:\n')
                    opName = t2[parse_int(opCode)]
                    dict(operations)[opName]()
    
                except KeyError, ValueError:
                    print KeyError
                    print ValueError
                    print 'Invalid Operation'
                except MenuReturn:
                    break
                
        except KeyboardInterrupt:
            sys.exit()
        except KeyError, ValueError:
            print 'Invalid Species Code!'
