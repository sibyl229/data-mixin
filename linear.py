from __future__ import division
import os
import numpy as np
from helpers import column_stack
from settings import *
#from DataRemap import remap

'''
Takes feature files and stack them along axis=1 (horizontally).
Note: same proteins and their order are expected across all feature files.
Finally, keep only proteins with values across all features. And normalize each feature the scores
'''

#
# take various feature files and vertically join them
#
pIDs = None
featuresStacked = None
hasValueStacked = None

for (dirpath, dirnames, filenames) in os.walk(FEATURE_FILE_PATH):
    for f in filenames:
        featureFilePath = os.path.join(dirpath, f)
        feature = np.genfromtxt(featureFilePath, delimiter='\t', 
            dtype=None, names=True)
        
        #
        # extract protein id, half-life, and features
        #
        idColName = 'Uniprot'
        excluded = ['Uniprot','hasValue'] # excluded from featuresStacked; other plans for these columns
        featureColNames = \
            [colName for colName in feature.dtype.names \
             if colName not in excluded]


        # stacking features
        if pIDs == None:
            pIDs = feature[['Uniprot']] # This is a structured array.
            featuresStacked = feature[featureColNames]
            allHasValue = feature['hasValue']
        else:
            if np.all(pIDs['Uniprot'] == feature['Uniprot']): 
                featuresStacked = column_stack(featuresStacked,
                                               feature[featureColNames])
                allHasValue = np.logical_and(allHasValue, feature['hasValue'])
            else:
                raise Exception('Mismatch protein ids. Please only combine features on the same list of proteins')


#
# take only protein entries, where all features are available
#
pIDs = pIDs[allHasValue]
featuresStacked = featuresStacked[allHasValue]

#
# backup the combined feature into single file
#
halfLifes = featuresStacked[['halflife_t12_in_h']]
otherColNames = \
    [cn for cn in featuresStacked.dtype.names if cn != 'halflife_t12_in_h']
combinedFeatures = column_stack(pIDs, halfLifes, 
                                featuresStacked[otherColNames])
combinedFeaturePath = os.path.join(CLEAN_INPUT_PATH,
                                   'combined_features.tsv')
with open(combinedFeaturePath, 'w') as f:
    np.savetxt(f, combinedFeatures, 
               delimiter='\t', comments='', fmt='%s',
               header='\t'.join(combinedFeatures.dtype.names))

#
# Normalized features, and backup
#
numOfFeatures = len(featuresStacked[0])
data = np.genfromtxt(combinedFeaturePath, 
                              dtype=float, skip_header=1,
                              usecols=np.arange(numOfFeatures)+1, #skip id column
                              delimiter='\t')
featuresArray = data[:,1:]
halfLifeArray = data[:,0]
means = np.mean(featuresArray, axis=0)
stds = np.std(featuresArray, axis=0)
normalizedFeatures = (featuresArray - means) / stds
normalizedFeaturePath = os.path.join(CLEAN_INPUT_PATH,
                                   'normalized_features.tsv')
headerNormalizedFeature = combinedFeatures.dtype.names[1:] # ignore id column for now
with open(normalizedFeaturePath, 'w') as f:
    np.savetxt(f, 
               np.hstack((halfLifeArray[:,None], normalizedFeatures)), 
               delimiter='\t', comments='', fmt='%s',
               header='\t'.join(headerNormalizedFeature)
    )


import pdb; pdb.set_trace()

