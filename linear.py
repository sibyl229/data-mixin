from __future__ import division
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
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
#ranges = np.max(featuresArray, axis=0) - np.min(featuresArray, axis=0)
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


#
# correlation between feature and HL
#
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)
subplots = [ax1, ax2, ax3, ax4]
for (j, featName) in enumerate(otherColNames):
    sbpl = subplots[j]
    # generate a heatmap (a 2d array) whose 1st dimension (vertical) represent different halfLife,
    # and 2nd dimension (horizontal) represents the feature score
    heatmap, xedges, yedges = \
        np.histogram2d(halfLifeArray, normalizedFeatures[:,j],
                       range=[[0,100],[-3,3]], 
                       bins=20)
    # due to skewness, normalized by total number of proteins of similar half-life (in the same bin)
    num_prot_by_halflife_bins = np.sum(heatmap,axis=1)[:,None] 
#    import pdb; pdb.set_trace()
    heatmap = heatmap / num_prot_by_halflife_bins
    heatmap = np.nan_to_num(heatmap)
    
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    # imshow display the 2darray like an image, so the heatmap needs to be rotated 
    # by 90 degree counterclockwise, for its first and second dimensions along x-axis and
    # y-axis respectively. 
    sbpl.imshow(np.rot90(heatmap), 
                interpolation="nearest", aspect='auto', extent=extent)
    sbpl.grid(True)
    sbpl.set_title(featName)
    sbpl.set_xlabel('Half-life in hours')
    sbpl.set_ylabel('Normalized feature score')

# sbpl.set_xlim(-1,1)
# sbpl.set_ylim(0,60)
#plt.show()


#
# fit linear model
#
ridge = Ridge()
ridge.fit(normalizedFeatures, halfLifeArray)
yHat = ridge.predict(normalizedFeatures)
print sp.stats.pearsonr(halfLifeArray, yHat)


fig, ax = plt.subplots()
ax.plot(halfLifeArray, yHat, 'o')
ax.set_ylim(0,60)
ax.set_xlim(0,60)
ax.set_title('Correlation between the real and the predicted')
#plt.show()

#import pdb; pdb.set_trace()

