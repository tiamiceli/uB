#!/usr/local/bin/python

#This code finds the BNB target center uncertainty from the pickle file from targetscan_ana.py
################################
# PREPARE ENVIRONMENT
################################
__author__ = 'miceli'
import sys
sys.path.insert(0, '/Users/tiamiceli/Documents/Code/PyBooNE')

import cPickle as pickle
import ROOT as rt
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib.colors import LogNorm
import myRootStyle
import scipy.stats as stats
from mpl_toolkits.mplot3d.axes3d import Axes3D
import pandas as pd
import csv

n_pulls=30
FinHFits=[]
FinHCovM=[]
SlugHFits=[]
SlugHCovM=[]
SlugVFits=[]
SlugVCovM=[]

print "loading file"
FinHFits = np.load("TargetFits_FinHFits.npy")
FinHCovM = np.load("TargetFits_FinHCovM.npy")
SlugHFits = np.load("TargetFits_SlugHFits.npy")
SlugHCovM = np.load("TargetFits_SlugHCovM.npy")
SlugVFits = np.load("TargetFits_SlugVFits.npy")
SlugVCovM = np.load("TargetFits_SlugVCovM.npy")
print "done loading"

TargetFits = [FinHFits, SlugHFits, SlugVFits]
TargetCovM = [FinHCovM, SlugHCovM, SlugVCovM]
bsname = ['Horizontal Fin Scan ', 'Horizontal Slug Scan ','Vertical Scan ']
pmname = ['BPM 875 ', 'BPM TG1 ', 'BPM TG2 ']
lmname = ['LM 875A ', 'LM 875B ', 'LM 875C ']

bs=0
for beamScan,beamScanCov in zip(TargetFits,TargetCovM):
    #print 'beamScan '+str(bs)
    pm=0
    for positionMon,positionMonCov in zip(beamScan,beamScanCov):
        #print '     positionMon '+str(pm)
        lm=0
        for lossMon,lossMonCov in zip(positionMon,positionMonCov):
            #print '         lossMon '+str(lm)

            leftPulls = np.random.multivariate_normal(lossMon[0],lossMonCov[0],n_pulls)
            midPulls = np.random.multivariate_normal(lossMon[1],lossMonCov[1],n_pulls)
            rightPulls = np.random.multivariate_normal(lossMon[2],lossMonCov[2],n_pulls)

            #with randompulls, iterate over them and find intersections with errors and central point with error
            #print '         find all roots'
            firstRoots = []
            for l in leftPulls:
                for m in midPulls:
                    root = (np.poly1d(m)-np.poly1d(l)).roots
                    firstRoots.append(root)

            secondRoots = []
            for m in midPulls:
                for r in rightPulls:
                    root = (np.poly1d(r)-np.poly1d(m)).roots
                    secondRoots.append(root)

            centers = []
            for r1 in firstRoots:
                for r2 in secondRoots:
                    center=np.mean([r1,r2])
                    centers.append(center)
            meanCenter = np.mean(centers)
            stdCenter = np.std(centers)

            print bsname[bs]+ pmname[pm]+lmname[lm]+ str(meanCenter) + " +/- " + str(stdCenter)
            #plt.hist(centers,50)
            #plt.show()

            lm=lm+1
        pm=pm+1
    bs=bs+1