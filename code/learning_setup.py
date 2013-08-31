from __future__ import division
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from helpers import column_stack
import go_helper
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
# keep only cytoplasmic proteins
#
cytoplasmic = go_helper.is_cytoplasmic(pIDs['Uniprot'])
pIDs = pIDs[cytoplasmic]
featuresStacked = featuresStacked[cytoplasmic]

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
maxes = np.max(featuresArray, axis=0)
mins = np.min(featuresArray, axis=0)
need_scaling = np.logical_and(maxes <= 1, mins >= -1)
scaledFeatures = featuresArray * 1
scaledFeatures[:,need_scaling] = \
    (featuresArray * 2)[:,need_scaling]

isbigValue = np.logical_or(featuresArray>3, featuresArray<-3)
bigValueRatio = np.sum(isbigValue, axis=0)/featuresArray.shape[0]
need_normalize = bigValueRatio > 0.05
# only normalize features with more than 5 percent fall out side of [-3,3]
normalizedFeatures = scaledFeatures * 1
normalizedFeatures[:,need_normalize] = \
    ((scaledFeatures - means) / stds)[:,need_normalize]
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
def featureVsHLHeatMap(halfLife, feature, normalizeBy='equi-HL-bin-size'):
    '''normalizeBy either 'equi-HL-bin-size' or 'equi-feature-bin-size' '''

    # generate a heatmap (a 2d array) whose 1st dimension (vertical) represent different halfLife,
    # and 2nd dimension (horizontal) represents the feature score
    heatmap, xedges, yedges = \
        np.histogram2d(halfLife, feature,
                       range=[[0,100],[-3,3]], 
                       bins=20)
    # due to skewness, normalized by total number of proteins of similar half-life (in the same bin)
    # or by total number of protein of similar feature value
    num_prot_by_halflife_bins = np.sum(heatmap,axis=1)[:,None] 
    num_prot_by_feature_bins = np.sum(heatmap,axis=0)[None,:] 
    if normalizeBy == 'equi-feature-bin-size':
        heatmap = heatmap / num_prot_by_feature_bins
    else:
        heatmap = heatmap / num_prot_by_halflife_bins
    heatmap = np.nan_to_num(heatmap)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap, extent


def plotFeatureAgainstHL(normalizedFeatureData, featureNames, normalizeBy='equi-feature-bin-size', numPerRow=2):
    numRows = int(np.ceil(len(featureNames) / numPerRow))
    fig, subplots = plt.subplots(numRows, numPerRow, sharex=True, sharey=True)

    #subplots = [ax1, ax2, ax3, ax4]
    for (j, featName) in enumerate(featureNames):
        if numRows > 1:
            sbpl = subplots[int(j/numPerRow)][ j % numPerRow]
        else:
            sbpl = subplots[j % numPerRow]

        heatmap, extent = featureVsHLHeatMap(halfLifeArray, 
                                             normalizedFeatureData[:,j], normalizeBy=normalizeBy)
        # imshow display the 2darray like an image, so the heatmap needs to be rotated 
        # by 90 degree counterclockwise, for its first and second dimensions along x-axis and
        # y-axis respectively. 
        sbpl.imshow(np.rot90(heatmap), 
                    interpolation="nearest", aspect='auto', extent=extent)
        sbpl.grid(True)
        sbpl.set_title(featName)
        #sbpl.set_xlabel('Half-life in hours')
        sbpl.set_ylabel('Norm. Score')#('Normalized feature score')

    return fig

def trivialCol(colArray):
    '''test if more than 99% of the column is 0'''
    numZeroes = np.sum(colArray==0)
    return numZeroes > len(colArray)*0.99

from matplotlib.backends.backend_pdf import PdfPages
if __name__ == '__main__':
    numRows = 2
    numPerRow = 2
    perSheet = numPerRow * numRows
    pp = PdfPages('fig/plot_features.pdf')
    pp2 = PdfPages('fig/plot_features2.pdf')
    #chosenColumns = np.array([not trivialCol(combinedFeatures[nm]) for nm in otherColNames])
    chosenColumns = np.array([True for nm in otherColNames])
    chosenColumnNames = np.array(otherColNames)[chosenColumns]
    chosenFeatures = normalizedFeatures[:,chosenColumns]
    #import pdb; pdb.set_trace()
    for i in range(0, len(chosenColumnNames), perSheet):
        header = chosenColumnNames[i : i+perSheet]
        data = chosenFeatures[:, i: i+perSheet]
        fig = plotFeatureAgainstHL(data, header, numPerRow=numPerRow, normalizeBy='equi-HL-bin-size')
        fig.savefig(pp, format='pdf')
        fig2 = plotFeatureAgainstHL(data, header, numPerRow=numPerRow, normalizeBy='equi-feature-bin-size')
        fig2.savefig(pp2, format='pdf')
        plt.close()
    pp.close()
    pp2.close()