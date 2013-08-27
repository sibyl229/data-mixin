from __future__ import division
from scipy import stats
import numpy as np
import sys
import pdb
import matplotlib.pyplot as plt

if len(sys.argv) >= 2:
	INPUT = sys.argv[1]
else: 
	INPUT = '/home/gw/sibyl/phalflife/input/15-60.tsv'
print INPUT
G = np.genfromtxt(INPUT, names=True)
numFeat = len(G.dtype.names)
c, r = 0, 0
P, R = np.zeros([numFeat,numFeat]), np.zeros([numFeat,numFeat])

for FEATURE in G.dtype.names:
	pcc = stats.pearsonr(G[FEATURE],G['halflife_t12_in_h'])
	print '#', FEATURE, pcc[0], pcc[1]
	r = 0
	for FEAT2 in G.dtype.names:
		pcc = stats.pearsonr(G[FEATURE],G[FEAT2])
		R[r,c] = pcc[0]
		P[r,c] = pcc[1]
		r = r+1
	c = c + 1

def makePlot(dataMatrix,OUTFILE):
	plt.clf()
	fig, ax = plt.subplots()
	image = dataMatrix
	ax.imshow(image, cmap=plt.cm.gray, interpolation='nearest')
	ax.set_title('PCC')
	# Move left and bottom spines outward by 10 points
	ax.spines['left'].set_position(('outward', 10))
	ax.spines['bottom'].set_position(('outward', 10))
	# Hide the right and top spines
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	# Only show ticks on the left and bottom spines
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.savefig(OUTFILE)
	return
makePlot(R, 'featCorrR.png')
makePlot(P,'featCorrP.png')
#makePlot(R*(P < 0.05), 'featCorrRCut.png')
