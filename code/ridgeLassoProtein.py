from __future__ import division
import numpy as np
from scipy import eye, stats
import matplotlib.pyplot as plt
from matplotlib import rc
from math import log, pow
import pylab
import pdb
#pdb.set_trace()
import argparse
import sys
import argparse

#parser = argparse.ArgumentParser(description='Command line options')



# read in data from default file or command line
if len(sys.argv) == 1: 
	sys.exit('ERROR: INPUT file missing, try: python ridgeLassoProtein.py INPUT')
else: 
	INPUT = sys.argv[1]

# remove N-end rule cols
if len(sys.argv) > 3 and int(sys.argv[3]) == 1:
	X = np.loadtxt(INPUT, skiprows=1, usecols=set(range(8)+range(28,34)) )
	G = np.genfromtxt(INPUT,usecols=set(range(8)+range(28,34)), names=True)
	thetaLabels = [head for head in G.dtype.names]
else:
	X = np.loadtxt(INPUT, skiprows=1)
	G = np.genfromtxt(INPUT, names=True)
	thetaLabels = [head for head in G.dtype.names]

thetaLabels = thetaLabels[1:]
np.random.shuffle(X)	# randomize data before splitting y and X
y = X[:,0] 	# response vector
X = X[:,1:]	# attributes: lcavol	lweight	age	lbph	svi	lcp	gleason	pgg45	lpsa
numFeat = np.shape(X)[1]

split = 0.5				# train size 50%, test size 50%
spRow = int(np.size(y)*split)
ytrain, ytest = y[0:spRow], y[spRow:] 	
Xtrain, Xtest = X[0:spRow], X[spRow:]

# compute standar score
# keep means and std of each attribute of X separately
Xbar = np.mean(Xtrain, axis=0) 
Xstd = np.std(Xtrain, axis=0) 
ybar = np.mean(ytrain) 
ytrain = ytrain - ybar
Xtrain = (Xtrain - Xbar) / Xstd
#Xtest = (Xtest - Xbar) / Xstd		# make (ytest,Xtest) standard score using mean and std of train
#ytest = ytest - ybar			# don't use train mean and std since if (ytest, Xtest) only 1 value makes no sense.

# calculate the solution to ridge: theta = (Xt*X + delta^2*I)^-1*Xt*y
# default value of delta is zero, least squares estimate


def ridge(X,y,d2=0):
	ridgeInv = np.linalg.inv(np.dot(X.T,X) + d2*eye(numFeat))
	theta = np.dot(ridgeInv.T,np.dot(X.T,y))
	return theta

d2Vals = [pow(10,myExp/10) for myExp in xrange(-20,40)]

# use labels for plot
# thetaLabels = ['f'+str(x+1) for x in range(numFeat)]

# question 1.1
if (1):
	d2BenchMark = [ridge(Xtrain,ytrain,d2) for d2 in d2Vals]
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xscale('log')
	plt.plot(d2Vals,d2BenchMark)
	plt.title('Ridge Regression Benchmark')
	plt.xlabel(r'$\delta^{2}$')
	plt.ylabel(r'$\theta$')
	plt.legend(thetaLabels, loc="center left")
	#ax.yaxis.label.set_size(40)
	#ax.xaxis.label.set_size(40)
	# colours http://stackoverflow.com/questions/8931268/using-colormaps-to-set-color-of-line-in-matplotlib
	#jet = cm = plt.get_cmap('jet') 
	#cNorm  = colors.Normalize(vmin=0, vmax=d2BenchMark)
	#scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	plt.savefig("ridgeThetaVsD2.png")
	#plt.show()

def relErr(yhat,y):
	diff = yhat - y
	relErr = np.dot(diff.T,diff)/np.dot(y.T,y)
	return relErr

def pearsonR(yhat,y):
	r = stats.pearsonr(yhat,y)[0]
	return r 

def yhat(Xstar, theta):		# Xstar should NOT be mean normalized, when use with train we need to undo mean normalization and get it back to original data values
	yhat = ybar + np.dot((Xstar - Xbar) / Xstd, theta)
	return yhat

def cv(d2Vals=[pow(10,myExp/10) for myExp in xrange(-20,40)],errFunc=relErr,trainMeth=ridge):	
		
	# compute relative error for test and train set: (y,x*) from these data sets
	testErrList = []
	trainErrList = []
	for d2 in d2Vals:
		theta = trainMeth(Xtrain,ytrain,d2)				# theta values always gotten the same way
		test_yhat = yhat(Xtest, theta)					# nx1 of predicted values for test set
		train_yhat = yhat(Xtrain*Xstd + Xbar, theta)			# have to get mean normalized Xtrain back to original format
	
		testErr = errFunc(test_yhat,ytest)			# compare to actual y values in test set and compute error
		trainErr = errFunc(train_yhat,ytrain + ybar)		# have to add back average that was previously subtracted off
		k=5
		diff = ytest - test_yhat
		testErrList.append(testErr)
		trainErrList.append(trainErr)
	return trainErrList, testErrList

# question 1.2
trainErrList, testErrList = cv(d2Vals, relErr,ridge)
if (1):
	# plot test and train error
	plotData=np.zeros(shape=(len(testErrList),2))
	plotData[:,0] = trainErrList
	plotData[:,1] = testErrList
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xscale('log')
	plt.plot(d2Vals,plotData)
	plt.title('Model Performance: Ridge')
	plt.xlabel(r'$\delta^{2}$')
	plt.ylabel(r'$\frac{\|\| y - X\theta \|\|_2^2}{\|\|y\|\|_2^2}$')
	plt.legend(['train','test'])
	plt.savefig("errVsD2Ridge.png")
	#plt.show()

# pearsonR
trainErrList, testErrList = cv(d2Vals, pearsonR,ridge)
if (1):
	# plot test and train error
	plotData=np.zeros(shape=(len(testErrList),2))
	plotData[:,0] = trainErrList
	plotData[:,1] = testErrList
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xscale('log')
	plt.plot(d2Vals,plotData)
	plt.title('Model Performance: Ridge')
	plt.xlabel(r'$\delta^{2}$')
	plt.ylabel(r'$PearsonCC$')
	plt.legend(['train','test'], loc="lower right")
	plt.savefig("rVsD2Ridge.png")
	#plt.show()

# question 1.3
def lasso(X, y, d2):
	row, col = X.shape
	d = col
	theta_ = np.zeros(col) 
	theta = ridge(X, y, d2) 
	c = np.zeros(d)
	a = 2 * np.dot(np.dot(X.T,X)*eye(d).T, np.ones(d))	# calculate aj outside cj loop
	while np.sum(np.abs(theta-theta_)) > 1e-5:
		theta_ = theta.copy()
		#c = [np.dot(2*Xtrain[:,j].T,ytrain - np.dot(Xtrain[:,:], theta) + Xtrain[:,j]*theta[j]) for j in range(d)]
		for j in range(d):
			c[j] = np.dot(2*Xtrain[:,j].T,ytrain - np.dot(Xtrain[:,:], theta) + Xtrain[:,j]*theta[j])
			if c[j] < -d2:
				theta[j] = (c[j] + d2) / a[j]
			elif c[j] > d2:
				theta[j] = (c[j] - d2) / a[j]
			else:
				theta[j] = 0	
	return theta

if (1):
	plotData = [lasso(Xtrain,ytrain,d2) for d2 in d2Vals]
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xscale('log')
	plt.plot(d2Vals,plotData)
	plt.title('Regularization Benchmark for Lasso')
	plt.legend(thetaLabels, loc="center left")
	plt.xlabel(r'$\delta^{2}$')
	plt.ylabel(r'$\theta$')
	plt.savefig("lassoThetaVsD2.png")
	#plt.show()


trainErrList, testErrList = cv(d2Vals, relErr,lasso)
if (1):
	# plot test and train error
	plotData=np.zeros(shape=(len(testErrList),2))
	plotData[:,0] = trainErrList
	plotData[:,1] = testErrList
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xscale('log')
	plt.plot(d2Vals,plotData)
	plt.title('Model Performance: Lasso ')
	plt.xlabel(r'$\delta^{2}$')
	plt.ylabel(r'$\frac{\|\| y - X\theta \|\|_2^2}{\|\|y\|\|_2^2}$')
	plt.legend(['train','test'])
	plt.savefig("errVsD2Lasso.png")
	#plt.show()

# pearsonR
trainErrList, testErrList = cv(d2Vals, pearsonR,lasso)
if (1):
	# plot test and train error
	plotData=np.zeros(shape=(len(testErrList),2))
	plotData[:,0] = trainErrList
	plotData[:,1] = testErrList
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xscale('log')
	plt.plot(d2Vals,plotData)
	plt.title('Model Performance: Lasso')
	plt.xlabel(r'$\delta^{2}$')
	plt.ylabel(r'$PearsonCC$')
	plt.legend(['train','test'], loc="lower right")
	plt.savefig("rVsD2Lasso.png")
	#plt.show()
if(1): 
	# choose d2 parameter, print out thetas for this parameter
	# choose theta
	if len(sys.argv) >= 3:
		d2 = float(sys.argv[2]) 
	else: d2 = 0
	theta = ridge(Xtrain,ytrain,d2)
	test_yhat = yhat(Xtest, theta)					# nx1 of predicted values for test set
	train_yhat = yhat(Xtrain*Xstd + Xbar, theta)			# have to get mean normalized Xtrain back to original format
	ytrain = ytrain + ybar
	for i,e, in enumerate(theta):
		if e != 0: print "#", thetaLabels[i], e
# OUTPUT
# lcavol 0.845534321278
# lweight 0.112863507676
# age -0.0953380205109
# lbph 0.067353621695
# svi 0.0669587916124
# gleason 0.0545111618812
#pdb.set_trace()
fig = plt.figure()
ax = plt.plot(ytrain,train_yhat,'x')
ax = plt.plot(ytest,test_yhat,'o')
plt.xlabel(r'$y$')
plt.ylabel(r'$yhat$')
plt.title('Model Predictions: Ridge for d2='+ str(d2))
plt.legend(['x: train','o: test'])
plt.savefig("predRidge.png")
