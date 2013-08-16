from __future__ import division
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from settings import *
import learning_setup


#
# fit linear model
#
normalizedFeatures = learning_setup.normalizedFeatures
halfLifeArray = learning_setup.halfLifeArray
ridge = Ridge()
ridge.fit(normalizedFeatures, halfLifeArray)
yHat = ridge.predict(normalizedFeatures)
print sp.stats.pearsonr(halfLifeArray, yHat)


ax = plt.subplot(aspect='equal')
ax.plot(halfLifeArray, yHat, 'o')
ax.set_ylim(0,60)
ax.set_xlim(0,60)
ax.set_title('Correlation between the real and the predicted')
plt.show()

#import pdb; pdb.set_trace()
'''
degr motif (Ken box)
number of Lysine
topfind (true n-end) or just remove M
take different length at N and C end (
yaakov levy weizmann
'''
