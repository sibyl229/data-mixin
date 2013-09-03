import sys
import os
from compute_data import *
from data import *
#from abc import ABCMeta, abstractmethod
from DataRemap import remap
from species.base import Species

    
def lookup_table(lst):
        return {i:sp for i,sp in enumerate(lst)}

def print_table(dtable, to_string=None):
    for key, item in dtable.items():
        if to_string:
            itemStr = to_string(item)
        else:
            itemStr = str(item)
        print '%3d\t%s' % (key, itemStr)

def parse_int(inStr):
    matchInt = re.search(r'(\d+)',inStr)
    intVal = None
    if matchInt:
        intVal = int(matchInt.group(1))
    return intVal
        

def compute_features(species):

    t = lookup_table(species.dataClasses)
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

def map_features(species):
    for dc in species.dataClasses:
        remap(dc())

#species_list = Species.ALL_SPECIES.values()


species_choices = [
    Species('species.human', 'Human'),
]

    
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
                    opDict = {
                        'Compute Features': (lambda: compute_features(species)),
                        'Map Features': (lambda: map_features(species)),
                        #   'Combine Features': lambda: ,
                    }
                    t2 = lookup_table(opDict.keys())
                    #import pdb; pdb.set_trace()

                    print_table(t2)
                    opCode = raw_input('Enter Operation Code:\n')
                    opName = t2[parse_int(opCode)]
                    opDict[opName]()
    
                except KeyError, ValueError:
                    print 'Invalid Operation'
                
        except KeyboardInterrupt:
            sys.exit()
        except KeyError, ValueError:
            print 'Invalid Species Code!'
