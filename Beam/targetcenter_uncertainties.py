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
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('TargetScanUncertainties_run0005544.pdf')

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

figFinH,  axarrFinH  = plt.subplots(3,3,figsize=(11,9),sharex=False,sharey=False)
figSlugH, axarrSlugH = plt.subplots(3,3,figsize=(11,9),sharex=False,sharey=False)
figSlugV, axarrSlugV = plt.subplots(3,3,figsize=(11,9),sharex=False,sharey=False)

fig = [figFinH, figSlugH, figSlugV]
axarr = [axarrFinH, axarrSlugH, axarrSlugV]

bs=0
for beamScan,beamScanCov in zip(TargetFits,TargetCovM):
    #print 'beamScan '+str(bs)
    fig[bs].subplots_adjust(left=0.1, bottom=0.1, top=0.92, wspace=0.3, hspace=0.3)

    fig[bs].text(0.5, 0.03, 'Center position (mm)', ha='center', size='large')
    fig[bs].text(0.015, 0.5, 'Frequency of center position', va='center', rotation='vertical', size='large')
    fig[bs].text(0.3,0.97, bsname[bs], size='large')
    axarr[bs][0,0].set_title(pmname[0])
    axarr[bs][0,1].set_title(pmname[1])
    axarr[bs][0,2].set_title(pmname[2])
    axarr[bs][0,0].set_ylabel(lmname[0])
    axarr[bs][1,0].set_ylabel(lmname[1])
    axarr[bs][2,0].set_ylabel(lmname[2])
    for l in range(0,3):
        for j in range(0,3):
            #f7axarr[f][l,j].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            axarr[bs][l,j].grid(True)



    pm=0
    for positionMon,positionMonCov in zip(beamScan,beamScanCov):
        #print '     positionMon '+str(pm)
        lm=0
        minx=999
        maxx=-999
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

            if np.amax(centers)>maxx:
                maxx=np.amax(centers)
            if np.amin(centers)>minx:
                minx=np.amin(centers)


            meanCenter = np.mean(centers)
            stdCenter = np.std(centers)

            string = bsname[bs]+ pmname[pm]+lmname[lm]+ str(meanCenter) + " +/- " + str(stdCenter)
            print string
            nums = str(meanCenter) + " +/- " + str(stdCenter)
            axarr[bs][lm,pm].hist(centers,50)
            #axarr[bs][lm,pm].annotate(nums)

            lm=lm+1

        #axarr[bs][0,pm].set_xlim(minx,maxx)
        #axarr[bs][1,pm].set_xlim(minx,maxx)
        #axarr[bs][2,pm].set_xlim(minx,maxx)

        pm=pm+1

    bs=bs+1


#plt.show()
pp.savefig(fig[0])
pp.savefig(fig[1])
pp.savefig(fig[2])
pp.close()