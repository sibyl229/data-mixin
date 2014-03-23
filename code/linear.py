from __future__ import division
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from settings import *
#import learning_setup


#
# fit linear model
#
def simple_linear(normalizedFeaturePath):
    data = np.genfromtxt(normalizedFeaturePath, dtype=None, skip_header=1, names=None, delimiter='\t')
    nonstatbleProteins = data[:,0] != 300
    data = data[nonstatbleProteins]
    halfLifes = data[:,0]
    normalizedFeatures = data[:,1:]

    ridge = Ridge()
    ridge.fit(normalizedFeatures, halfLifes)
    yHat = ridge.predict(normalizedFeatures)
    print data.shape
    print sp.stats.pearsonr(halfLifes, yHat)

    ax = plt.subplot(aspect='equal')
    ax.plot(halfLifes, yHat, 'o')
    ax.set_ylim(0,200)
    ax.set_xlim(0,200)
    ax.set_xlabel('Actual Half-life')
    ax.set_ylabel('Predicted Half-life')
    ax.set_title('Correlation between the real and the predicted')
    plt.show()


