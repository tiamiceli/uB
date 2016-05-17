#!/usr/local/bin/python

#This code finds the BNB target center based on horizontal and vertical scans

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
pp = PdfPages('TargetScan_run0005544.pdf')

#Force this style on all histograms
mysty = myRootStyle.myRootStyle()
rt.gROOT.SetStyle("mysty")
rt.gROOT.ForceStyle();

################################
# READ IN LOCATION INFORMATION
################################
inputloc = open('/Users/tiamiceli/Data/Beam/locationsBeamInstrumentation.csv','rb')
reader = csv.reader(inputloc)
# eliminate blank rows if they exist
rows = [row for row in reader if row]
location = {}
for (key, value) in rows:
    value = float(value)
    location.setdefault(key,[value])

#print "\n location['HT860D'] \n"
#print location['HT860D'][0]

################################
# READ IN ROOT FILE TREE
################################

tree = rt.TChain("bnb")
tree.Add("/Users/tiamiceli/Data/Beam/beamdata_0005544_000.root")
tree.Add("/Users/tiamiceli/Data/Beam/beamdata_0005544_001.root")
treeEntries = tree.GetEntries()
nonzeroEntries = 0
################################
# ARRAYS TO MAKE PLOTS
################################

atimestamp=np.zeros(treeEntries)

aTOR875=np.zeros(treeEntries)

aHI860=np.zeros(treeEntries)
aHI861=np.zeros(treeEntries)
aHI864=np.zeros(treeEntries)
aHI866=np.zeros(treeEntries)
aHI868=np.zeros(treeEntries)
aHI870=np.zeros(treeEntries)
aHI872=np.zeros(treeEntries)
aHI873=np.zeros(treeEntries)
aHI875=np.zeros(treeEntries)
aHITG1=np.zeros(treeEntries)
aHITG2=np.zeros(treeEntries)

aHP860=np.zeros(treeEntries)
aHP861=np.zeros(treeEntries)
aHP864=np.zeros(treeEntries)
aHP866=np.zeros(treeEntries)
aHP868=np.zeros(treeEntries)
aHP870=np.zeros(treeEntries)
aHP872=np.zeros(treeEntries)
aHP873=np.zeros(treeEntries)
aHP875=np.zeros(treeEntries)
aHPTG1=np.zeros(treeEntries)
aHPTG2=np.zeros(treeEntries)

aVI860=np.zeros(treeEntries)
aVI861=np.zeros(treeEntries)
aVI864=np.zeros(treeEntries)
aVI867=np.zeros(treeEntries)
aVI869=np.zeros(treeEntries)
aVI870=np.zeros(treeEntries)
aVI871=np.zeros(treeEntries)
aVI873=np.zeros(treeEntries)
aVI875=np.zeros(treeEntries)
aVITG1=np.zeros(treeEntries)
aVITG2=np.zeros(treeEntries)
    
aVP860=np.zeros(treeEntries)
aVP861=np.zeros(treeEntries)
aVP864=np.zeros(treeEntries)
aVP867=np.zeros(treeEntries)
aVP869=np.zeros(treeEntries)
aVP870=np.zeros(treeEntries)
aVP871=np.zeros(treeEntries)
aVP873=np.zeros(treeEntries)
aVP875=np.zeros(treeEntries)
aVPTG1=np.zeros(treeEntries)
aVPTG2=np.zeros(treeEntries)

aM875BB_mean_x=np.zeros(treeEntries)
aM875BB_rms_x =np.zeros(treeEntries)
aM875BB_mean_y=np.zeros(treeEntries)
aM875BB_rms_y =np.zeros(treeEntries)
aM876BB_mean_x=np.zeros(treeEntries)
aM876BB_rms_x =np.zeros(treeEntries)
aM876BB_mean_y=np.zeros(treeEntries)
aM876BB_rms_y =np.zeros(treeEntries)
aMMBTBB_mean_x=np.zeros(treeEntries)
aMMBTBB_rms_x =np.zeros(treeEntries)
aMMBTBB_mean_y=np.zeros(treeEntries)
aMMBTBB_rms_y =np.zeros(treeEntries)

aLM851A=np.zeros(treeEntries)
aLM860 =np.zeros(treeEntries)
aLM862 =np.zeros(treeEntries)
aLM864 =np.zeros(treeEntries)
aLM865A=np.zeros(treeEntries)
aLM865B=np.zeros(treeEntries)
aLM866 =np.zeros(treeEntries)
aLM867 =np.zeros(treeEntries)
aLM868 =np.zeros(treeEntries)
aLM869 =np.zeros(treeEntries)
aLM870 =np.zeros(treeEntries)
aLM871 =np.zeros(treeEntries)
aLM872 =np.zeros(treeEntries)
aLM873 =np.zeros(treeEntries)
aLM874 =np.zeros(treeEntries)
aLM875A=np.zeros(treeEntries)
aLM875B=np.zeros(treeEntries)
aLM875C=np.zeros(treeEntries)
aLMBPST=np.zeros(treeEntries)

aLM875ATOR=np.zeros(treeEntries)
aLM875BTOR=np.zeros(treeEntries)
aLM875CTOR=np.zeros(treeEntries)

scanType=np.zeros(treeEntries)

#scanIdx = [np.array(),np.array(),np.array(),np.array()]

################################
# LOOP OVER EVENTS
################################

for ientry in xrange(treeEntries):
    nb = tree.GetEntry(ientry)
    if nb <= 0:
        continue
    if (ientry%500==0):
        print "ientry = "+str(ientry)+"/"+str(treeEntries)

    missingInfo = tree.HI875<-0.9 or tree.HITG1<-0.9 or tree.HITG2<-0.9 or\
                       tree.VI875<-0.9 or tree.VITG1<-0.9 or tree.VITG2<-0.9 or\
                       tree.HP875<-200 or tree.HPTG1<-200 or tree.HPTG2<-200 or\
                       tree.VP875<-200 or tree.VPTG1<-200 or tree.VPTG2<-200 or\
                       tree.VI873<-0.9 or tree.VI864<-0.9 or tree.HI873<-0.9 or\
                       tree.VI860<-0.9 or tree.HI870<-0.9 or tree.HI864<-0.9 or\
                       tree.VP864<-200 or tree.VP869<-200 or tree.VP870<-200 or tree.TOR875<4.3


        #continue
    #else, fill up my numpy arrays!

    time = tree.timestamp - 1458671333542
    atimestamp[nonzeroEntries]=time
    if time > 150000 and time < 375000:
        scanType[nonzeroEntries] = 1 #Vertical Scan
    #    scanIdx[1].append(nonzeroEntries)
    elif time > 375000 and time < 460000:
        scanType[nonzeroEntries] = 2 #Horizontal Fin Scan
    #    scanIdx[2].append(nonzeroEntries)
    elif time > 500000 and time < 725000:
        scanType[nonzeroEntries] = 3 #Horizontal Slug Scan
    #    scanIdx[3].append(nonzeroEntries)
    #else:
    #    scanIdx[0].append(nonzeroEntries)
    if missingInfo:
        scanType[nonzeroEntries]=-1

    aTOR875[nonzeroEntries]=tree.TOR875

    aHI860[nonzeroEntries]=tree.HI860
    aHI861[nonzeroEntries]=tree.HI861
    aHI864[nonzeroEntries]=tree.HI864
    aHI866[nonzeroEntries]=tree.HI866
    aHI868[nonzeroEntries]=tree.HI868
    aHI870[nonzeroEntries]=tree.HI870
    aHI872[nonzeroEntries]=tree.HI872
    aHI873[nonzeroEntries]=tree.HI873
    aHI875[nonzeroEntries]=tree.HI875
    aHITG1[nonzeroEntries]=tree.HITG1
    aHITG2[nonzeroEntries]=tree.HITG2
    
    aHP860[nonzeroEntries]=tree.HP860
    aHP861[nonzeroEntries]=tree.HP861
    aHP864[nonzeroEntries]=tree.HP864
    aHP866[nonzeroEntries]=tree.HP866
    aHP868[nonzeroEntries]=tree.HP868
    aHP870[nonzeroEntries]=tree.HP870
    aHP872[nonzeroEntries]=tree.HP872
    aHP873[nonzeroEntries]=tree.HP873
    aHP875[nonzeroEntries]=-(tree.HP875)
    aHPTG1[nonzeroEntries]=tree.HPTG1
    aHPTG2[nonzeroEntries]=tree.HPTG2

    aVI860[nonzeroEntries]=tree.VI860
    aVI861[nonzeroEntries]=tree.VI861
    aVI864[nonzeroEntries]=tree.VI864
    aVI867[nonzeroEntries]=tree.VI867
    aVI869[nonzeroEntries]=tree.VI869
    aVI870[nonzeroEntries]=tree.VI870
    aVI871[nonzeroEntries]=tree.VI871
    aVI873[nonzeroEntries]=tree.VI873
    aVI875[nonzeroEntries]=tree.VI875
    aVITG1[nonzeroEntries]=tree.VITG1
    aVITG2[nonzeroEntries]=tree.VITG2
    
    aVP860[nonzeroEntries]=tree.VP860
    aVP861[nonzeroEntries]=tree.VP861
    aVP864[nonzeroEntries]=tree.VP864
    aVP867[nonzeroEntries]=tree.VP867
    aVP869[nonzeroEntries]=tree.VP869
    aVP870[nonzeroEntries]=tree.VP870
    aVP871[nonzeroEntries]=tree.VP871
    aVP873[nonzeroEntries]=tree.VP873
    aVP875[nonzeroEntries]=tree.VP875
    aVPTG1[nonzeroEntries]=tree.VPTG1
    aVPTG2[nonzeroEntries]=tree.VPTG2

    aM875BB_mean_x[nonzeroEntries]=tree.M875BB_mean_x
    aM875BB_rms_x[nonzeroEntries] =tree.M875BB_rms_x
    aM875BB_mean_y[nonzeroEntries]=tree.M875BB_mean_y
    aM875BB_rms_y[nonzeroEntries] =tree.M875BB_rms_y
    aM876BB_mean_x[nonzeroEntries]=tree.M876BB_mean_x
    aM876BB_rms_x[nonzeroEntries] =tree.M876BB_rms_x
    aM876BB_mean_y[nonzeroEntries]=tree.M876BB_mean_y
    aM876BB_rms_y[nonzeroEntries] =tree.M876BB_rms_y
    aMMBTBB_mean_x[nonzeroEntries]=tree.MMBTBB_mean_x
    aMMBTBB_rms_x[nonzeroEntries] =tree.MMBTBB_rms_x
    aMMBTBB_mean_y[nonzeroEntries]=tree.MMBTBB_mean_y
    aMMBTBB_rms_y[nonzeroEntries] =tree.MMBTBB_rms_y

    aLM851A[nonzeroEntries]=tree.LM851A
    aLM860[nonzeroEntries] =tree.LM860
    aLM862[nonzeroEntries] =tree.LM862
    aLM864[nonzeroEntries] =tree.LM864
    aLM865A[nonzeroEntries]=tree.LM865A
    aLM865B[nonzeroEntries]=tree.LM865B
    aLM866[nonzeroEntries] =tree.LM866
    aLM867[nonzeroEntries] =tree.LM867
    aLM868[nonzeroEntries] =tree.LM868
    aLM869[nonzeroEntries] =tree.LM869
    aLM870[nonzeroEntries] =tree.LM870
    aLM871[nonzeroEntries] =tree.LM871
    aLM872[nonzeroEntries] =tree.LM872
    aLM873[nonzeroEntries] =tree.LM873
    aLM874[nonzeroEntries] =tree.LM874
    aLM875A[nonzeroEntries]=tree.LM875A
    aLM875B[nonzeroEntries]=tree.LM875B
    aLM875C[nonzeroEntries]=tree.LM875C
    aLMBPST[nonzeroEntries]=tree.LMBPST

    aLM875ATOR[nonzeroEntries]= tree.LM875A/tree.TOR875
    aLM875BTOR[nonzeroEntries]= tree.LM875B/tree.TOR875
    aLM875CTOR[nonzeroEntries]= tree.LM875C/tree.TOR875

    nonzeroEntries = nonzeroEntries + 1

#end for loop over events


for i in range(-1,4):
    print 'scanType=='+str(i)+' is size: ' + str(len(scanType[np.where(scanType==i)]))
name = ['Regular Run ', 'Vertical Scan ', 'Horizontal Fin Scan ', 'Horizontal Slug Scan ']
################################
# MAKE PLOTS
################################
FinHFits =[]
FinHCovM =[]
SlugHFits=[]
SlugHCovM=[]
SlugVFits=[]
SlugVCovM=[]

#---------------plot position vs. time---------------
fig2, f2axarr = plt.subplots(2,2,figsize=(11,9))
fig2.subplots_adjust(wspace=0.3, hspace=0.3)
f2axarr[0,0].scatter(atimestamp[np.where(scanType!=-1)], aHPTG2[np.where(scanType!=-1)], c='teal'  , alpha=0.1, lw=0, label='HPTG2')
f2axarr[0,0].set_xlim([np.amin(atimestamp[np.where(scanType!=-1)]),np.amax(atimestamp[np.where(scanType!=-1)])]);
f2axarr[0,0].set_ylim([np.amin(aHPTG2[np.where(scanType!=-1)]),np.amax(aHPTG2[np.where(scanType!=-1)])]);
f2axarr[0,0].set_ylabel('Horizontal Position (mm)')
f2axarr[0,0].set_xlabel('Time (a.u.)')
f2axarr[0,0].ticklabel_format(style='sci', axis='x', scilimits=(0,5))
f2axarr[0,0].minorticks_on()
f2axarr[0,0].grid(True)
f2axarr[0,0].add_patch(pat.Rectangle([150000,-4.5],225000,12,fill=False,ec='red'))
f2axarr[0,0].add_patch(pat.Rectangle([375000,-4.5],85000,12,fill=False,ec='red'))
f2axarr[0,0].add_patch(pat.Rectangle([500000,-4.5],225000,12,fill=False,ec='red'))

f2axarr[1,0].scatter(atimestamp[np.where(scanType!=-1)], aVPTG2[np.where(scanType!=-1)], c='violet'  , alpha=0.1, lw=0, label='VPTG2')
f2axarr[1,0].set_xlim([np.amin(atimestamp[np.where(scanType!=-1)]),np.amax(atimestamp[np.where(scanType!=-1)])]);
f2axarr[1,0].set_ylim([np.amin(aVPTG2[np.where(scanType!=-1)]),np.amax(aVPTG2[np.where(scanType!=-1)])]);
f2axarr[1,0].set_ylabel('Vertical Position (mm)')
f2axarr[1,0].set_xlabel('Time (a.u.)')
f2axarr[1,0].ticklabel_format(style='sci', axis='x', scilimits=(0,5))
f2axarr[1,0].minorticks_on()
f2axarr[1,0].grid(True)
f2axarr[1,0].add_patch(pat.Rectangle([150000,-9.5],225000,19.5,fill=False,ec='red'))
f2axarr[1,0].add_patch(pat.Rectangle([375000,-9.5],85000,19.5,fill=False,ec='red'))
f2axarr[1,0].add_patch(pat.Rectangle([500000,-9.5],225000,19.5,fill=False,ec='red'))

f2axarr[0,1].scatter(aHPTG2[np.where(scanType!=-1)], aVPTG2[np.where(scanType!=-1)], c='orange'  , alpha=0.1, lw=0, label='VPTG2')
f2axarr[0,1].set_xlim([np.amin(aHPTG2[np.where(scanType!=-1)]),np.amax(aHPTG2[np.where(scanType!=-1)])]);
f2axarr[0,1].set_ylim([np.amin(aVPTG2[np.where(scanType!=-1)]),np.amax(aVPTG2[np.where(scanType!=-1)])]);
f2axarr[0,1].set_xlabel('Horizontal Position (mm)')
f2axarr[0,1].set_ylabel('Vertical Position (mm)')
f2axarr[0,1].minorticks_on()
f2axarr[0,1].grid(True)

f2axarr[1,1].scatter(atimestamp[np.where(scanType!=-1)], aTOR875[np.where(scanType!=-1)], c='steelblue'  , alpha=0.8, lw=0, label='TOR875')
#f2axarr[1,1].set_xlim([np.amin(aHPTG2[np.where(scanType!=-1)]),np.amax(aHPTG2[np.where(scanType!=-1)])]);
#f2axarr[1,1].set_ylim(2.5,5);
f2axarr[1,1].set_xlabel('Time (a.u.)')
f2axarr[1,1].set_ylabel('Toroid 875')
f2axarr[1,1].ticklabel_format(style='sci', axis='x', scilimits=(0,5))
f2axarr[1,1].minorticks_on()
f2axarr[1,1].grid(True)
f2axarr[1,1].add_patch(pat.Rectangle([150000,4.47],225000,0.2,fill=False,ec='red'))
f2axarr[1,1].add_patch(pat.Rectangle([375000,4.47],85000,0.2,fill=False,ec='red'))
f2axarr[1,1].add_patch(pat.Rectangle([500000,4.47],225000,0.2,fill=False,ec='red'))


fig2.text(0.2,0.97, 'Positions and Times at HPTG2 and VPTG2 & no missing info', size='large')

#---------------plot vertical position vs. beam line position---------------
fig0_0, a0_0 = plt.subplots(1,2,figsize=(10,6))
fig0_1, a0_1 = plt.subplots(1,2,figsize=(10,6))
fig0_2, a0_2 = plt.subplots(1,2,figsize=(10,6))
#fig0_3, a0_3 = plt.subplots(1,2,figsize=(10,6))
fig0 = [fig0_0, fig0_1, fig0_2]#, fig0_3]
a0 = [a0_0, a0_1, a0_2]#, a0_3]

f=0
for i in range(1,4):
    fig0[f].subplots_adjust(right=0.8)
    a0[f][0].scatter(aVP860[np.where(scanType==i)], location['VP860'][0] * np.ones(len(aVP860[np.where(scanType==i)])), c='firebrick', alpha=0.4, lw=0, label='VP860')
    a0[f][0].scatter(aVP861[np.where(scanType==i)], location['VP861'][0] * np.ones(len(aVP861[np.where(scanType==i)])), c='red', alpha=0.4, lw=0, label='VP861')
    a0[f][0].scatter(aVP864[np.where(scanType==i)], location['VP864'][0] * np.ones(len(aVP864[np.where(scanType==i)])), c='orange', alpha=0.4, lw=0, label='VP864')
    a0[f][0].scatter(aVP867[np.where(scanType==i)], location['VP867'][0] * np.ones(len(aVP867[np.where(scanType==i)])), c='gold', alpha=0.4, lw=0, label='VP867')
    a0[f][0].scatter(aVP869[np.where(scanType==i)], location['VP869'][0] * np.ones(len(aVP869[np.where(scanType==i)])), c='lawngreen', alpha=0.4, lw=0, label='VP869')
    a0[f][0].scatter(aVP870[np.where(scanType==i)], location['VP870'][0] * np.ones(len(aVP870[np.where(scanType==i)])), c='green', alpha=0.4, lw=0, label='VP870')
    a0[f][0].scatter(aVP871[np.where(scanType==i)], location['VP871'][0] * np.ones(len(aVP871[np.where(scanType==i)])), c='teal', alpha=0.4, lw=0, label='VP871')
    a0[f][0].scatter(aVP873[np.where(scanType==i)], location['VP873'][0] * np.ones(len(aVP873[np.where(scanType==i)])), c='blue', alpha=0.4, lw=0, label='VP873')
    a0[f][0].scatter(aVP875[np.where(scanType==i)], location['VP875'][0] * np.ones(len(aVP875[np.where(scanType==i)])), c='violet', alpha=0.4, lw=0, label='VP875')
    a0[f][0].scatter(aVPTG1[np.where(scanType==i)], location['VPTG1'][0] * np.ones(len(aVPTG1[np.where(scanType==i)])), c='purple', alpha=0.4, lw=0, label='VPTG1')
    a0[f][0].scatter(aVPTG2[np.where(scanType==i)], location['VPTG2'][0] * np.ones(len(aVPTG2[np.where(scanType==i)])), c='magenta', alpha=0.4, lw=0, label='VPTG2')
    a0[f][0].minorticks_on()
    a0[f][0].grid(True)

    a0[f][0].set_xlabel('Vertical Position (mm)')
    a0[f][0].set_ylabel('Longitudinal Position (m)')
    fig0[f].text(0.3,0.97, name[i]+'Beam Position v. Longitudinal Position', size='large')
    #a0[f].ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    a0[f][1].scatter((aHP860[np.where(scanType==i)]), location['HP860'][0] * np.ones(len(aHP860[np.where(scanType==i)])), c='firebrick', alpha=0.4, lw=0, label='HP860')
    a0[f][1].scatter(aHP861[np.where(scanType==i)], location['HP861'][0] * np.ones(len(aHP861[np.where(scanType==i)])), c='red', alpha=0.4, lw=0, label='HP861')
    a0[f][1].scatter(aHP864[np.where(scanType==i)], location['HP864'][0] * np.ones(len(aHP864[np.where(scanType==i)])), c='orange', alpha=0.4, lw=0, label='HP864')
    a0[f][1].scatter(aHP866[np.where(scanType==i)], location['HP866'][0] * np.ones(len(aHP866[np.where(scanType==i)])), c='gold', alpha=0.4, lw=0, label='HP866')
    a0[f][1].scatter(aHP868[np.where(scanType==i)], location['HP868'][0] * np.ones(len(aHP868[np.where(scanType==i)])), c='lawngreen', alpha=0.4, lw=0, label='HP868')
    a0[f][1].scatter(aHP870[np.where(scanType==i)], location['HP870'][0] * np.ones(len(aHP870[np.where(scanType==i)])), c='green', alpha=0.4, lw=0, label='HP870')
    a0[f][1].scatter(aHP872[np.where(scanType==i)], location['HP872'][0] * np.ones(len(aHP872[np.where(scanType==i)])), c='teal', alpha=0.4, lw=0, label='HP872')
    a0[f][1].scatter(aHP873[np.where(scanType==i)], location['HP873'][0] * np.ones(len(aHP873[np.where(scanType==i)])), c='blue', alpha=0.4, lw=0, label='HP873')
    a0[f][1].scatter(aHP875[np.where(scanType==i)], location['HP875'][0] * np.ones(len(aHP875[np.where(scanType==i)])), c='violet', alpha=0.4, lw=0, label='HP875\'')
    a0[f][1].scatter(aHPTG1[np.where(scanType==i)], location['HPTG1'][0] * np.ones(len(aHPTG1[np.where(scanType==i)])), c='purple', alpha=0.4, lw=0, label='HPTG1')
    a0[f][1].scatter(aHPTG2[np.where(scanType==i)], location['HPTG2'][0] * np.ones(len(aHPTG2[np.where(scanType==i)])), c='magenta', alpha=0.4, lw=0, label='HPTG2')
    a0[f][1].minorticks_on()
    a0[f][1].grid(True)

    a0[f][1].set_xlabel('Horizontal Position (mm)')
    a0[f][1].set_ylabel('Longitudinal Position (m)')
    a0[f][1].legend(loc='center right' ,bbox_to_anchor=(0.98,0.5), bbox_transform=fig0[f].transFigure, scatterpoints=1)
    f=f+1

#---------------plot Loss vs. horizontal beam position ---------------
fig6_1, f6axarr_1 = plt.subplots(3,4,figsize=(11,9),sharex=True,sharey=True)
fig6_2, f6axarr_2 = plt.subplots(3,4,figsize=(11,9),sharex=True,sharey=True)
#fig6_3, f6axarr_3 = plt.subplots(3,4,figsize=(11,9),sharex=True,sharey=True)
#fig6_4, f6axarr_4 = plt.subplots(3,4,figsize=(11,9),sharex=True,sharey=True)

fig6 = [fig6_1, fig6_2]#, fig6_3, fig6_4]
f6axarr = [f6axarr_1, f6axarr_2]#, f6axarr_3, f6axarr_4]
#plt.subplots_adjust(right=0.8)
#plt.minorticks_on()
f=0
for i in range(2,4): ##for i in range(0,4):
    f6axarr[f][0,3].set_ylim(0,0.04)
    f6axarr[f][0,3].set_xlim(-8,12)
    f6axarr[f][0,3].scatter(aHP875[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875A')
    f6axarr[f][0,3].scatter(aHPTG1[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875A')
    f6axarr[f][0,3].scatter(aHPTG2[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875A')
    f6axarr[f][0,0].scatter(aHP875[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875A')
    f6axarr[f][0,1].scatter(aHPTG1[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875A')
    f6axarr[f][0,2].scatter(aHPTG2[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875A')

    f6axarr[f][1,3].scatter(aHP875[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875B')
    f6axarr[f][1,3].scatter(aHPTG1[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875B')
    f6axarr[f][1,3].scatter(aHPTG2[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875B')
    f6axarr[f][1,0].scatter(aHP875[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875B')
    f6axarr[f][1,1].scatter(aHPTG1[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875B')
    f6axarr[f][1,2].scatter(aHPTG2[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875B')

    f6axarr[f][2,3].scatter(aHP875[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875C')
    f6axarr[f][2,3].scatter(aHPTG1[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875C')
    f6axarr[f][2,3].scatter(aHPTG2[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875C')
    f6axarr[f][2,0].scatter(aHP875[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875C')
    f6axarr[f][2,1].scatter(aHPTG1[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875C')
    f6axarr[f][2,2].scatter(aHPTG2[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875C')

    f6axarr[f][0,0].set_ylabel('LM875A')
    f6axarr[f][1,0].set_ylabel('LM875B')
    f6axarr[f][2,0].set_ylabel('LM875C')
    f6axarr[f][0,0].set_title('HP875\'')
    f6axarr[f][0,1].set_title('HPTG1')
    f6axarr[f][0,2].set_title('HPTG2')
    f6axarr[f][0,3].set_title('overlay')

#    f6axarr[f][2,3].legend(loc='center right',bbox_to_anchor=(0.99,0.5),bbox_transform=plt.gcf().transFigure, scatterpoints=1)

    #fig6.tight_layout()
    fig6[f].subplots_adjust(left=0.1, bottom=0.1, top=0.92, wspace=0.3, hspace=0.3)

    fig6[f].text(0.5, 0.04, 'Position (Horizontal) (mm)', ha='center', size='large')
    fig6[f].text(0.015, 0.5, 'Loss/TOR875 (a.u.)', va='center', rotation='vertical', size='large')
    fig6[f].text(0.3,0.97, name[i]+'Loss/TOR875 v. Horizontal Position', size='large')

    for l in range(0,3):
        for j in range(0,4):

            #f6axarr[f][l,j].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            f6axarr[f][l,j].grid(True)
            
    ################
    #  FIT  LINES  #
    ################
    #i==2 #fin scan
    firstCol= [2,3,   3,4,   4,7]
    firstColA=[0.5,2.5]
    #firstCol= [2,3,   3,4.25,   4.5,7]
    #firstColA=[0.5,2.5]
    secondCol=[-1.25,-0.25, -0.25,0.5,     0.5,3  ]
    secondColA=[-2,-0.5]
    thirdCol= [-1,0,  0,0.75,     0.75,3  ]
    thirdColA= [-1.5,0]
    if i==3: #slug scan
        firstCol= [-2.25,-0.5,   0,7,   8,10]
        firstColA=[-4,-2]
        secondCol=[-4,-3,     -2,3,  4,5 ]
        secondColA=[secondCol[0],secondCol[1]]
        thirdCol= [-3.5,-2,     -2,3,  4,5 ]
        thirdColA= [-4.5,-3]



    # HP875
    xHP875l  = aHP875[  np.where((scanType==i) & (aHP875>firstCol[0]) & (aHP875<firstCol[1]))]
    xHP875lA = aHP875[  np.where((scanType==i) & (aHP875>firstColA[0]) & (aHP875<firstColA[1]))]
    y1l875 = aLM875ATOR[np.where((scanType==i) & (aHP875>firstColA[0]) & (aHP875<firstColA[1]))]
    y2l875 = aLM875BTOR[np.where((scanType==i) & (aHP875>firstCol[0]) & (aHP875<firstCol[1]))]
    y3l875 = aLM875CTOR[np.where((scanType==i) & (aHP875>firstCol[0]) & (aHP875<firstCol[1]))]
    yl875 = [y1l875, y2l875, y3l875]
    xHP875t  = aHP875[  np.where((scanType==i) & (aHP875>firstCol[2]) & (aHP875<firstCol[3]))]
    y1t875 = aLM875ATOR[np.where((scanType==i) & (aHP875>firstCol[2]) & (aHP875<firstCol[3]))]
    y2t875 = aLM875BTOR[np.where((scanType==i) & (aHP875>firstCol[2]) & (aHP875<firstCol[3]))]
    y3t875 = aLM875CTOR[np.where((scanType==i) & (aHP875>firstCol[2]) & (aHP875<firstCol[3]))]
    yt875 = [y1t875, y2t875, y3t875]
    xHP875r  = aHP875[  np.where((scanType==i) & (aHP875>firstCol[4]) & (aHP875<firstCol[5]))]
    y1r875 = aLM875ATOR[np.where((scanType==i) & (aHP875>firstCol[4]) & (aHP875<firstCol[5]))]
    y2r875 = aLM875BTOR[np.where((scanType==i) & (aHP875>firstCol[4]) & (aHP875<firstCol[5]))]
    y3r875 = aLM875CTOR[np.where((scanType==i) & (aHP875>firstCol[4]) & (aHP875<firstCol[5]))]
    yr875 = [y1r875, y2r875, y3r875]
    xHP875 = [xHP875l, xHP875t, xHP875r]
    y875 = [yl875, yt875, yr875]
    
    fitvectorA875 = [np.poly1d(np.polyfit(xHP875lA ,y875[0][0],deg=1)), np.poly1d(np.polyfit(xHP875[1],y875[1][0],deg=1)), np.poly1d(np.polyfit(xHP875[2],y875[2][0],deg=1))]
    fitvectorB875 = [np.poly1d(np.polyfit(xHP875[0],y875[0][1],deg=1)), np.poly1d(np.polyfit(xHP875[1],y875[1][1],deg=1)), np.poly1d(np.polyfit(xHP875[2],y875[2][1],deg=1))]
    fitvectorC875 = [np.poly1d(np.polyfit(xHP875[0],y875[0][2],deg=1)), np.poly1d(np.polyfit(xHP875[1],y875[1][2],deg=1)), np.poly1d(np.polyfit(xHP875[2],y875[2][2],deg=1))]
    fitvector875 = [fitvectorA875, fitvectorB875, fitvectorC875]

    polyfitA875L, polyCovMA875L =np.polyfit(xHP875lA ,y875[0][0],deg=1,cov=True)
    polyfitA875M, polyCovMA875M =np.polyfit(xHP875[1],y875[1][0],deg=1,cov=True)
    polyfitA875R, polyCovMA875R=np.polyfit(xHP875[2],y875[2][0],deg=1,cov=True)
    polyfitA875=[polyfitA875L,polyfitA875M,polyfitA875R]
    polyCovMA875=[polyCovMA875L,polyCovMA875M,polyCovMA875R]

    polyfitB875L, polyCovMB875L =np.polyfit(xHP875[0],y875[0][1],deg=1,cov=True)
    polyfitB875M, polyCovMB875M =np.polyfit(xHP875[1],y875[1][1],deg=1,cov=True)
    polyfitB875R, polyCovMB875R =np.polyfit(xHP875[2],y875[2][1],deg=1,cov=True)
    polyfitB875 = [polyfitB875L,polyfitB875M,polyfitB875R]
    polyCovMB875= [polyCovMB875L,polyCovMB875M,polyCovMB875R]

    polyfitC875L, polyCovMC875L =np.polyfit(xHP875[0],y875[0][2],deg=1,cov=True)
    polyfitC875M, polyCovMC875M =np.polyfit(xHP875[1],y875[1][2],deg=1,cov=True)
    polyfitC875R, polyCovMC875R =np.polyfit(xHP875[2],y875[2][2],deg=1,cov=True)
    polyfitC875=[polyfitC875L,polyfitC875M,polyfitC875R]
    polyCovMC875 = [polyCovMC875L,polyCovMC875M,polyCovMC875R]
    
    polyfitvector875 = [polyfitA875,polyfitB875,polyfitC875]
    polyfitCovM875= [polyCovMA875,polyCovMB875,polyCovMC875]

    #if i==2 #fin scan
    linspace875 = [np.linspace(-7,5,6), np.linspace(1,5,6), np.linspace(1,10,6)]
    if i==3: #slug scan
        linspace875 = [np.linspace(-7,2,6), np.linspace(-5,10,6), np.linspace(6,10,6)]

    # HPTG1
    xHPTG1l  = aHPTG1[  np.where((scanType==i) & (aHPTG1>secondCol[0]) & (aHPTG1<secondCol[1]))]
    xHPTG1lA  = aHPTG1[  np.where((scanType==i) & (aHPTG1>secondColA[0]) & (aHPTG1<secondColA[1]))]
    y1lTG1 = aLM875ATOR[np.where((scanType==i) & (aHPTG1>secondColA[0]) & (aHPTG1<secondColA[1]))]
    y2lTG1 = aLM875BTOR[np.where((scanType==i) & (aHPTG1>secondCol[0]) & (aHPTG1<secondCol[1]))]
    y3lTG1 = aLM875CTOR[np.where((scanType==i) & (aHPTG1>secondCol[0]) & (aHPTG1<secondCol[1]))]
    ylTG1 = [y1lTG1, y2lTG1, y3lTG1]
    xHPTG1t  = aHPTG1[  np.where((scanType==i) & (aHPTG1>secondCol[2]) & (aHPTG1<secondCol[3]))]
    y1tTG1 = aLM875ATOR[np.where((scanType==i) & (aHPTG1>secondCol[2]) & (aHPTG1<secondCol[3]))]
    y2tTG1 = aLM875BTOR[np.where((scanType==i) & (aHPTG1>secondCol[2]) & (aHPTG1<secondCol[3]))]
    y3tTG1 = aLM875CTOR[np.where((scanType==i) & (aHPTG1>secondCol[2]) & (aHPTG1<secondCol[3]))]
    ytTG1 = [y1tTG1, y2tTG1, y3tTG1]
    xHPTG1r  = aHPTG1[  np.where((scanType==i) & (aHPTG1>secondCol[4]) & (aHPTG1<secondCol[5]))]
    y1rTG1 = aLM875ATOR[np.where((scanType==i) & (aHPTG1>secondCol[4]) & (aHPTG1<secondCol[5]))]
    y2rTG1 = aLM875BTOR[np.where((scanType==i) & (aHPTG1>secondCol[4]) & (aHPTG1<secondCol[5]))]
    y3rTG1 = aLM875CTOR[np.where((scanType==i) & (aHPTG1>secondCol[4]) & (aHPTG1<secondCol[5]))]
    yrTG1 = [y1rTG1, y2rTG1, y3rTG1]
    xHPTG1 = [xHPTG1l, xHPTG1t, xHPTG1r]
    yTG1 = [ylTG1, ytTG1, yrTG1]

    fitvectorATG1 = [np.poly1d(np.polyfit(xHPTG1lA ,yTG1[0][0],deg=1)), np.poly1d(np.polyfit(xHPTG1[1],yTG1[1][0],deg=1)), np.poly1d(np.polyfit(xHPTG1[2],yTG1[2][0],deg=1))]
    fitvectorBTG1 = [np.poly1d(np.polyfit(xHPTG1[0],yTG1[0][1],deg=1)), np.poly1d(np.polyfit(xHPTG1[1],yTG1[1][1],deg=1)), np.poly1d(np.polyfit(xHPTG1[2],yTG1[2][1],deg=1))]
    fitvectorCTG1 = [np.poly1d(np.polyfit(xHPTG1[0],yTG1[0][2],deg=1)), np.poly1d(np.polyfit(xHPTG1[1],yTG1[1][2],deg=1)), np.poly1d(np.polyfit(xHPTG1[2],yTG1[2][2],deg=1))]
    fitvectorTG1 = [fitvectorATG1, fitvectorBTG1, fitvectorCTG1]

    
    polyfitATG1L, polyCovMATG1L =np.polyfit(xHPTG1lA ,yTG1[0][0],deg=1,cov=True)
    polyfitATG1M, polyCovMATG1M =np.polyfit(xHPTG1[1],yTG1[1][0],deg=1,cov=True)
    polyfitATG1R, polyCovMATG1R=np.polyfit(xHPTG1[2],yTG1[2][0],deg=1,cov=True)
    polyfitATG1=[polyfitATG1L,polyfitATG1M,polyfitATG1R]
    polyCovMATG1=[polyCovMATG1L,polyCovMATG1M,polyCovMATG1R]

    polyfitBTG1L, polyCovMBTG1L =np.polyfit(xHPTG1[0],yTG1[0][1],deg=1,cov=True)
    polyfitBTG1M, polyCovMBTG1M =np.polyfit(xHPTG1[1],yTG1[1][1],deg=1,cov=True)
    polyfitBTG1R, polyCovMBTG1R =np.polyfit(xHPTG1[2],yTG1[2][1],deg=1,cov=True)
    polyfitBTG1 = [polyfitBTG1L,polyfitBTG1M,polyfitBTG1R]
    polyCovMBTG1= [polyCovMBTG1L,polyCovMBTG1M,polyCovMBTG1R]

    polyfitCTG1L, polyCovMCTG1L =np.polyfit(xHPTG1[0],yTG1[0][2],deg=1,cov=True)
    polyfitCTG1M, polyCovMCTG1M =np.polyfit(xHPTG1[1],yTG1[1][2],deg=1,cov=True)
    polyfitCTG1R, polyCovMCTG1R =np.polyfit(xHPTG1[2],yTG1[2][2],deg=1,cov=True)
    polyfitCTG1=[polyfitCTG1L,polyfitCTG1M,polyfitCTG1R]
    polyCovMCTG1 = [polyCovMCTG1L,polyCovMCTG1M,polyCovMCTG1R]
    
    polyfitvectorTG1 = [polyfitATG1,polyfitBTG1,polyfitCTG1]
    polyfitCovMTG1= [polyCovMATG1,polyCovMBTG1,polyCovMCTG1]
    

    #if i==2 #fin scan
    linspaceTG1 = [np.linspace(-7,2,6), np.linspace(-2,2,6), np.linspace(-1,7,6)]
    if i==3: #slug scan
        linspaceTG1 = [np.linspace(-7,-1,6), np.linspace(-5,5.5,6), np.linspace(2,7,6)]

    # HPTG2
    xHPTG2l  = aHPTG2[  np.where((scanType==i) & (aHPTG2>thirdCol[0]) & (aHPTG2<thirdCol[1]))]
    xHPTG2lA  = aHPTG2[  np.where((scanType==i) & (aHPTG2>thirdColA[0]) & (aHPTG2<thirdColA[1]))]
    y1lTG2 = aLM875ATOR[np.where((scanType==i) & (aHPTG2>thirdColA[0]) & (aHPTG2<thirdColA[1]))]
    y2lTG2 = aLM875BTOR[np.where((scanType==i) & (aHPTG2>thirdCol[0]) & (aHPTG2<thirdCol[1]))]
    y3lTG2 = aLM875CTOR[np.where((scanType==i) & (aHPTG2>thirdCol[0]) & (aHPTG2<thirdCol[1]))]
    ylTG2 = [y1lTG2, y2lTG2, y3lTG2]
    xHPTG2t  = aHPTG2[  np.where((scanType==i) & (aHPTG2>thirdCol[2]) & (aHPTG2<thirdCol[3]))]
    y1tTG2 = aLM875ATOR[np.where((scanType==i) & (aHPTG2>thirdCol[2]) & (aHPTG2<thirdCol[3]))]
    y2tTG2 = aLM875BTOR[np.where((scanType==i) & (aHPTG2>thirdCol[2]) & (aHPTG2<thirdCol[3]))]
    y3tTG2 = aLM875CTOR[np.where((scanType==i) & (aHPTG2>thirdCol[2]) & (aHPTG2<thirdCol[3]))]
    ytTG2 = [y1tTG2, y2tTG2, y3tTG2]
    xHPTG2r  = aHPTG2[  np.where((scanType==i) & (aHPTG2>thirdCol[4]) & (aHPTG2<thirdCol[5]))]
    y1rTG2 = aLM875ATOR[np.where((scanType==i) & (aHPTG2>thirdCol[4]) & (aHPTG2<thirdCol[5]))]
    y2rTG2 = aLM875BTOR[np.where((scanType==i) & (aHPTG2>thirdCol[4]) & (aHPTG2<thirdCol[5]))]
    y3rTG2 = aLM875CTOR[np.where((scanType==i) & (aHPTG2>thirdCol[4]) & (aHPTG2<thirdCol[5]))]
    yrTG2 = [y1rTG2, y2rTG2, y3rTG2]
    xHPTG2 = [xHPTG2l, xHPTG2t, xHPTG2r]
    yTG2 = [ylTG2, ytTG2, yrTG2]

    fitvectorATG2 = [np.poly1d(np.polyfit(xHPTG2lA,yTG2[0][0],deg=1)), np.poly1d(np.polyfit(xHPTG2[1],yTG2[1][0],deg=1)), np.poly1d(np.polyfit(xHPTG2[2],yTG2[2][0],deg=1))]
    fitvectorBTG2 = [np.poly1d(np.polyfit(xHPTG2[0],yTG2[0][1],deg=1)), np.poly1d(np.polyfit(xHPTG2[1],yTG2[1][1],deg=1)), np.poly1d(np.polyfit(xHPTG2[2],yTG2[2][1],deg=1))]
    fitvectorCTG2 = [np.poly1d(np.polyfit(xHPTG2[0],yTG2[0][2],deg=1)), np.poly1d(np.polyfit(xHPTG2[1],yTG2[1][2],deg=1)), np.poly1d(np.polyfit(xHPTG2[2],yTG2[2][2],deg=1))]
    fitvectorTG2 = [fitvectorATG2, fitvectorBTG2, fitvectorCTG2]

    polyfitATG2L, polyCovMATG2L =np.polyfit(xHPTG2lA ,yTG2[0][0],deg=1,cov=True)
    #print "********************************************"
    #print "********************************************"
    #print "********************************************"
    #print "xHPTG2lA.shape"
    #print xHPTG2lA.shape
    #print "yTG2[0][0].shape"
    #print yTG2[0][0].shape
    #print "polyfitATG2L.shape"
    #print polyfitATG2L.shape
    #print "polyfitATG2L"
    #print polyfitATG2L
    #print "********************************************"
    #print "********************************************"
    #print "********************************************"
    polyfitATG2M, polyCovMATG2M =np.polyfit(xHPTG2[1],yTG2[1][0],deg=1,cov=True)
    polyfitATG2R, polyCovMATG2R=np.polyfit(xHPTG2[2],yTG2[2][0],deg=1,cov=True)
    polyfitATG2=[polyfitATG2L,polyfitATG2M,polyfitATG2R]
    polyCovMATG2=[polyCovMATG2L,polyCovMATG2M,polyCovMATG2R]

    polyfitBTG2L, polyCovMBTG2L =np.polyfit(xHPTG2[0],yTG2[0][1],deg=1,cov=True)
    polyfitBTG2M, polyCovMBTG2M =np.polyfit(xHPTG2[1],yTG2[1][1],deg=1,cov=True)
    polyfitBTG2R, polyCovMBTG2R =np.polyfit(xHPTG2[2],yTG2[2][1],deg=1,cov=True)
    polyfitBTG2 = [polyfitBTG2L,polyfitBTG2M,polyfitBTG2R]
    polyCovMBTG2= [polyCovMBTG2L,polyCovMBTG2M,polyCovMBTG2R]

    polyfitCTG2L, polyCovMCTG2L =np.polyfit(xHPTG2[0],yTG2[0][2],deg=1,cov=True)
    polyfitCTG2M, polyCovMCTG2M =np.polyfit(xHPTG2[1],yTG2[1][2],deg=1,cov=True)
    polyfitCTG2R, polyCovMCTG2R =np.polyfit(xHPTG2[2],yTG2[2][2],deg=1,cov=True)
    polyfitCTG2=[polyfitCTG2L,polyfitCTG2M,polyfitCTG2R]
    polyCovMCTG2 = [polyCovMCTG2L,polyCovMCTG2M,polyCovMCTG2R]

    
    polyfitvectorTG2 = [polyfitATG2,polyfitBTG2,polyfitCTG2]
    polyfitCovMTG2= [polyCovMATG2,polyCovMBTG2,polyCovMCTG2]

    #if i==2 #fin scan
    linspaceTG2 = [np.linspace(-7,3,6), np.linspace(-1,3,6), np.linspace(-1,7,6)]
    if i==3: #slug scan
        linspaceTG2 = [np.linspace(-7,0,6), np.linspace(-5,5.5,6), np.linspace(2,7,6)]

    fitvector = [fitvector875, fitvectorTG1, fitvectorTG2]
    linspace = [linspace875, linspaceTG1, linspaceTG2]
    #print '******HORIZONAL SCAN*******'

    for pm, posMon in enumerate(fitvector):
        for lm,lossMon in enumerate(posMon):
            #print 'plot['+str(lm)+','+str(pm)+'] '+ name[i]
            roots = []
            for sec,section in enumerate(linspace[pm]):
                f6axarr[f][lm,pm].plot(section,lossMon[sec](section),'k-',alpha=0.5)
                if sec>0:
                    root = (lossMon[sec]-lossMon[sec-1]).roots
                    f6axarr[f][lm,pm].plot(root[0],lossMon[sec](root[0]),'o')
                    roots.append(root[0])

            #print 'roots = '
            #print roots

            center = float(roots[0]+(0.5*(roots[1]-roots[0])))
            centery = float(lossMon[1](center))
            f6axarr[f][lm,pm].plot(center,centery,'o')
            annotation_string = r"%0.2f" % (center)
            f6axarr[f][lm,pm].annotate(annotation_string, (center-2, centery-0.005), color='red')

    if i==2: #fin
        FinHFits = [polyfitvector875, polyfitvectorTG1, polyfitvectorTG2]
        FinHCovM = [polyfitCovM875, polyfitCovMTG1, polyfitCovMTG2]
    
    if i==3: #slug
        SlugHFits = [polyfitvector875, polyfitvectorTG1, polyfitvectorTG2]
        SlugHCovM = [polyfitCovM875, polyfitCovMTG1, polyfitCovMTG2]
        
    
    f=f+1

#---------------plot Loss vs. vertical beam position ---------------
fig7_1, f7axarr_1 = plt.subplots(3,4,figsize=(11,9),sharex=True,sharey=True)
#fig7_2, f7axarr_2 = plt.subplots(3,4,figsize=(11,9),sharex=True,sharey=True)
#fig7_3, f7axarr_3 = plt.subplots(3,4,figsize=(11,9),sharex=True,sharey=True)
#fig7_4, f7axarr_4 = plt.subplots(3,4,figsize=(11,9),sharex=True,sharey=True)

fig7 = [fig7_1]#, fig7_2, fig7_3, fig7_4]
f7axarr = [f7axarr_1]#, f7axarr_2, f7axarr_3, f7axarr_4]

#plt.subplots_adjust(right=0.8)
#plt.minorticks_on()

f=0
for i in range(1,2): ##for i in range(0,4):
    f7axarr[f][0,3].set_ylim(0,0.04)
    f7axarr[f][0,3].set_xlim(-8,7)
    f7axarr[f][0,3].scatter(aVP875[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875A')
    f7axarr[f][0,3].scatter(aVPTG1[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875A')
    f7axarr[f][0,3].scatter(aVPTG2[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875A')
    f7axarr[f][0,0].scatter(aVP875[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875A')
    f7axarr[f][0,1].scatter(aVPTG1[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875A')
    f7axarr[f][0,2].scatter(aVPTG2[np.where(scanType==i)], aLM875ATOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875A')

    f7axarr[f][1,3].scatter(aVP875[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875B')
    f7axarr[f][1,3].scatter(aVPTG1[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875B')
    f7axarr[f][1,3].scatter(aVPTG2[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875B')
    f7axarr[f][1,0].scatter(aVP875[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875B')
    f7axarr[f][1,1].scatter(aVPTG1[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875B')
    f7axarr[f][1,2].scatter(aVPTG2[np.where(scanType==i)], aLM875BTOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875B')

    f7axarr[f][2,3].scatter(aVP875[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875C')
    f7axarr[f][2,3].scatter(aVPTG1[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875C')
    f7axarr[f][2,3].scatter(aVPTG2[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875C')
    f7axarr[f][2,0].scatter(aVP875[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='firebrick'  , alpha=0.1, lw=0, label='LM875C')
    f7axarr[f][2,1].scatter(aVPTG1[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='orange'     , alpha=0.1, lw=0, label='LM875C')
    f7axarr[f][2,2].scatter(aVPTG2[np.where(scanType==i)], aLM875CTOR[np.where(scanType==i)], c='steelblue'  , alpha=0.1, lw=0, label='LM875C')



    f7axarr[f][0,0].set_ylabel('LM875A')
    f7axarr[f][1,0].set_ylabel('LM875B')
    f7axarr[f][2,0].set_ylabel('LM875C')
    f7axarr[f][0,0].set_title('VP875')
    f7axarr[f][0,1].set_title('VPTG1')
    f7axarr[f][0,2].set_title('VPTG2')
    f7axarr[f][0,3].set_title('overlay')

    #f7axarr[f][2,3].legend(loc='center right',bbox_to_anchor=(0.99,0.5),bbox_transform=plt.gcf().transFigure, scatterpoints=1)

    #fig7[f].tight_layout()
    fig7[f].subplots_adjust(left=0.1, bottom=0.1, top=0.92, wspace=0.3, hspace=0.3)

    fig7[f].text(0.5, 0.03, 'Position (Vertical) (mm)', ha='center', size='large')
    fig7[f].text(0.015, 0.5, 'Loss/TOR875 (a.u.)', va='center', rotation='vertical', size='large')
    fig7[f].text(0.3,0.97, name[i]+'Loss/TOR875 v. Vertical Position', size='large')

    for l in range(0,3):
        for j in range(0,4):
            #f7axarr[f][l,j].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            f7axarr[f][l,j].grid(True)

    ################
    #  FIT  LINES  #
    ################

    # VP875
    xVP875l  = aVP875[    np.where((scanType==i) & (aVP875>-7) & (aVP875<-4))]
    y1l875 = aLM875ATOR[np.where((scanType==i) & (aVP875>-7) & (aVP875<-4))]
    y2l875 = aLM875BTOR[np.where((scanType==i) & (aVP875>-7) & (aVP875<-4))]
    y3l875 = aLM875CTOR[np.where((scanType==i) & (aVP875>-7) & (aVP875<-4))]
    yl875 = [y1l875, y2l875, y3l875]
    xVP875t  = aVP875[    np.where((scanType==i) & (aVP875>-3) & (aVP875<3))]
    y1t875 = aLM875ATOR[np.where((scanType==i) & (aVP875>-3) & (aVP875<3))]
    y2t875 = aLM875BTOR[np.where((scanType==i) & (aVP875>-3) & (aVP875<3))]
    y3t875 = aLM875CTOR[np.where((scanType==i) & (aVP875>-3) & (aVP875<3))]
    yt875 = [y1t875, y2t875, y3t875]
    xVP875r  = aVP875[    np.where((scanType==i) & (aVP875>4) & (aVP875<7))]
    y1r875 = aLM875ATOR[np.where((scanType==i) & (aVP875>4) & (aVP875<7))]
    y2r875 = aLM875BTOR[np.where((scanType==i) & (aVP875>4) & (aVP875<7))]
    y3r875 = aLM875CTOR[np.where((scanType==i) & (aVP875>4) & (aVP875<7))]
    yr875 = [y1r875, y2r875, y3r875]
    xVP875 = [xVP875l, xVP875t, xVP875r]
    y875 = [yl875, yt875, yr875]
    
    fitvectorA875 = [np.poly1d(np.polyfit(xVP875[0],y875[0][0],deg=1)), np.poly1d(np.polyfit(xVP875[1],y875[1][0],deg=1)), np.poly1d(np.polyfit(xVP875[2],y875[2][0],deg=1))]
    fitvectorB875 = [np.poly1d(np.polyfit(xVP875[0],y875[0][1],deg=1)), np.poly1d(np.polyfit(xVP875[1],y875[1][1],deg=1)), np.poly1d(np.polyfit(xVP875[2],y875[2][1],deg=1))]
    fitvectorC875 = [np.poly1d(np.polyfit(xVP875[0],y875[0][2],deg=1)), np.poly1d(np.polyfit(xVP875[1],y875[1][2],deg=1)), np.poly1d(np.polyfit(xVP875[2],y875[2][2],deg=1))]
    fitvector875 = [fitvectorA875, fitvectorB875, fitvectorC875]
    linspace875 = [np.linspace(-7,-2,6), np.linspace(-5,5,6), np.linspace(2,7,6)]

    polyfitA875L, polyCovMA875L =np.polyfit(xVP875[0] ,y875[0][0],deg=1,cov=True)
    polyfitA875M, polyCovMA875M =np.polyfit(xVP875[1],y875[1][0],deg=1,cov=True)
    polyfitA875R, polyCovMA875R= np.polyfit(xVP875[2],y875[2][0],deg=1,cov=True)
    polyfitA875=[polyfitA875L,polyfitA875M,polyfitA875R]
    polyCovMA875=[polyCovMA875L,polyCovMA875M,polyCovMA875R]

    polyfitB875L, polyCovMB875L =np.polyfit(xVP875[0],y875[0][1],deg=1,cov=True)
    polyfitB875M, polyCovMB875M =np.polyfit(xVP875[1],y875[1][1],deg=1,cov=True)
    polyfitB875R, polyCovMB875R =np.polyfit(xVP875[2],y875[2][1],deg=1,cov=True)
    polyfitB875 = [polyfitB875L,polyfitB875M,polyfitB875R]
    polyCovMB875= [polyCovMB875L,polyCovMB875M,polyCovMB875R]

    polyfitC875L, polyCovMC875L =np.polyfit(xVP875[0],y875[0][2],deg=1,cov=True)
    polyfitC875M, polyCovMC875M =np.polyfit(xVP875[1],y875[1][2],deg=1,cov=True)
    polyfitC875R, polyCovMC875R =np.polyfit(xVP875[2],y875[2][2],deg=1,cov=True)
    polyfitC875=[polyfitC875L,polyfitC875M,polyfitC875R]
    polyCovMC875 = [polyCovMC875L,polyCovMC875M,polyCovMC875R]
    
    polyfitvector875 = [polyfitA875,polyfitB875,polyfitC875]
    polyfitCovM875= [polyCovMA875,polyCovMB875,polyCovMC875]

    # VPTG1
    xVPTG1l  = aVPTG1[    np.where((scanType==i) & (aVPTG1>-7) & (aVPTG1<-4.25))]
    y1lTG1 = aLM875ATOR[np.where((scanType==i) & (aVPTG1>-7) & (aVPTG1<-4.25))]
    y2lTG1 = aLM875BTOR[np.where((scanType==i) & (aVPTG1>-7) & (aVPTG1<-4.25))]
    y3lTG1 = aLM875CTOR[np.where((scanType==i) & (aVPTG1>-7) & (aVPTG1<-4.25))]
    ylTG1 = [y1lTG1, y2lTG1, y3lTG1]
    xVPTG1t  = aVPTG1[    np.where((scanType==i) & (aVPTG1>-3) & (aVPTG1<1))]
    y1tTG1 = aLM875ATOR[np.where((scanType==i) & (aVPTG1>-3) & (aVPTG1<1))]
    y2tTG1 = aLM875BTOR[np.where((scanType==i) & (aVPTG1>-3) & (aVPTG1<1))]
    y3tTG1 = aLM875CTOR[np.where((scanType==i) & (aVPTG1>-3) & (aVPTG1<1))]
    ytTG1 = [y1tTG1, y2tTG1, y3tTG1]
    xVPTG1r  = aVPTG1[    np.where((scanType==i) & (aVPTG1>2.5) & (aVPTG1<4.5))]
    y1rTG1 = aLM875ATOR[np.where((scanType==i) & (aVPTG1>2.5) & (aVPTG1<4.5))]
    y2rTG1 = aLM875BTOR[np.where((scanType==i) & (aVPTG1>2.5) & (aVPTG1<4.5))]
    y3rTG1 = aLM875CTOR[np.where((scanType==i) & (aVPTG1>2.5) & (aVPTG1<4.5))]
    yrTG1 = [y1rTG1, y2rTG1, y3rTG1]
    xVPTG1 = [xVPTG1l, xVPTG1t, xVPTG1r]
    yTG1 = [ylTG1, ytTG1, yrTG1]

    fitvectorATG1 = [np.poly1d(np.polyfit(xVPTG1[0],yTG1[0][0],deg=1)), np.poly1d(np.polyfit(xVPTG1[1],yTG1[1][0],deg=1)), np.poly1d(np.polyfit(xVPTG1[2],yTG1[2][0],deg=1))]
    fitvectorBTG1 = [np.poly1d(np.polyfit(xVPTG1[0],yTG1[0][1],deg=1)), np.poly1d(np.polyfit(xVPTG1[1],yTG1[1][1],deg=1)), np.poly1d(np.polyfit(xVPTG1[2],yTG1[2][1],deg=1))]
    fitvectorCTG1 = [np.poly1d(np.polyfit(xVPTG1[0],yTG1[0][2],deg=1)), np.poly1d(np.polyfit(xVPTG1[1],yTG1[1][2],deg=1)), np.poly1d(np.polyfit(xVPTG1[2],yTG1[2][2],deg=1))]
    fitvectorTG1 = [fitvectorATG1, fitvectorBTG1, fitvectorCTG1]
    linspaceTG1 = [np.linspace(-7,-2,6), np.linspace(-5,5,6), np.linspace(0,7,6)]

    polyfitATG1L, polyCovMATG1L =np.polyfit(xVPTG1[0] ,yTG1[0][0],deg=1,cov=True)
    polyfitATG1M, polyCovMATG1M =np.polyfit(xVPTG1[1],yTG1[1][0],deg=1,cov=True)
    polyfitATG1R, polyCovMATG1R= np.polyfit(xVPTG1[2],yTG1[2][0],deg=1,cov=True)
    polyfitATG1=[polyfitATG1L,polyfitATG1M,polyfitATG1R]
    polyCovMATG1=[polyCovMATG1L,polyCovMATG1M,polyCovMATG1R]

    polyfitBTG1L, polyCovMBTG1L =np.polyfit(xVPTG1[0],yTG1[0][1],deg=1,cov=True)
    polyfitBTG1M, polyCovMBTG1M =np.polyfit(xVPTG1[1],yTG1[1][1],deg=1,cov=True)
    polyfitBTG1R, polyCovMBTG1R =np.polyfit(xVPTG1[2],yTG1[2][1],deg=1,cov=True)
    polyfitBTG1 = [polyfitBTG1L,polyfitBTG1M,polyfitBTG1R]
    polyCovMBTG1= [polyCovMBTG1L,polyCovMBTG1M,polyCovMBTG1R]

    polyfitCTG1L, polyCovMCTG1L =np.polyfit(xVPTG1[0],yTG1[0][2],deg=1,cov=True)
    polyfitCTG1M, polyCovMCTG1M =np.polyfit(xVPTG1[1],yTG1[1][2],deg=1,cov=True)
    polyfitCTG1R, polyCovMCTG1R =np.polyfit(xVPTG1[2],yTG1[2][2],deg=1,cov=True)
    polyfitCTG1=[polyfitCTG1L,polyfitCTG1M,polyfitCTG1R]
    polyCovMCTG1 = [polyCovMCTG1L,polyCovMCTG1M,polyCovMCTG1R]
   
    polyfitvectorTG1 = [polyfitATG1,polyfitBTG1,polyfitCTG1]
    polyfitCovMTG1= [polyCovMATG1,polyCovMBTG1,polyCovMCTG1]


    # VPTG2
    xVPTG2l  = aVPTG2[    np.where((scanType==i) & (aVPTG2>-6) & (aVPTG2<-4))]
    y1lTG2 = aLM875ATOR[np.where((scanType==i) & (aVPTG2>-6) & (aVPTG2<-4))]
    y2lTG2 = aLM875BTOR[np.where((scanType==i) & (aVPTG2>-6) & (aVPTG2<-4))]
    y3lTG2 = aLM875CTOR[np.where((scanType==i) & (aVPTG2>-6) & (aVPTG2<-4))]
    ylTG2 = [y1lTG2, y2lTG2, y3lTG2]
    xVPTG2t  = aVPTG2[    np.where((scanType==i) & (aVPTG2>-2) & (aVPTG2<2))]
    y1tTG2 = aLM875ATOR[np.where((scanType==i) & (aVPTG2>-2) & (aVPTG2<2))]
    y2tTG2 = aLM875BTOR[np.where((scanType==i) & (aVPTG2>-2) & (aVPTG2<2))]
    y3tTG2 = aLM875CTOR[np.where((scanType==i) & (aVPTG2>-2) & (aVPTG2<2))]
    ytTG2 = [y1tTG2, y2tTG2, y3tTG2]
    xVPTG2r  = aVPTG2[    np.where((scanType==i) & (aVPTG2>3) & (aVPTG2<5))]
    y1rTG2 = aLM875ATOR[np.where((scanType==i) & (aVPTG2>3) & (aVPTG2<5))]
    y2rTG2 = aLM875BTOR[np.where((scanType==i) & (aVPTG2>3) & (aVPTG2<5))]
    y3rTG2 = aLM875CTOR[np.where((scanType==i) & (aVPTG2>3) & (aVPTG2<5))]
    yrTG2 = [y1rTG2, y2rTG2, y3rTG2]
    xVPTG2 = [xVPTG2l, xVPTG2t, xVPTG2r]
    yTG2 = [ylTG2, ytTG2, yrTG2]

    fitvectorATG2 = [np.poly1d(np.polyfit(xVPTG2[0],yTG2[0][0],deg=1)), np.poly1d(np.polyfit(xVPTG2[1],yTG2[1][0],deg=1)), np.poly1d(np.polyfit(xVPTG2[2],yTG2[2][0],deg=1))]
    fitvectorBTG2 = [np.poly1d(np.polyfit(xVPTG2[0],yTG2[0][1],deg=1)), np.poly1d(np.polyfit(xVPTG2[1],yTG2[1][1],deg=1)), np.poly1d(np.polyfit(xVPTG2[2],yTG2[2][1],deg=1))]
    fitvectorCTG2 = [np.poly1d(np.polyfit(xVPTG2[0],yTG2[0][2],deg=1)), np.poly1d(np.polyfit(xVPTG2[1],yTG2[1][2],deg=1)), np.poly1d(np.polyfit(xVPTG2[2],yTG2[2][2],deg=1))]
    fitvectorTG2 = [fitvectorATG2, fitvectorBTG2, fitvectorCTG2]
    linspaceTG2 = [np.linspace(-7,-2,6), np.linspace(-5,5,6), np.linspace(2,7,6)]
    
    polyfitATG2L, polyCovMATG2L =np.polyfit(xVPTG2[0] ,yTG2[0][0],deg=1,cov=True)
    polyfitATG2M, polyCovMATG2M =np.polyfit(xVPTG2[1],yTG2[1][0],deg=1,cov=True)
    polyfitATG2R, polyCovMATG2R= np.polyfit(xVPTG2[2],yTG2[2][0],deg=1,cov=True)
    polyfitATG2=[polyfitATG2L,polyfitATG2M,polyfitATG2R]
    polyCovMATG2=[polyCovMATG2L,polyCovMATG2M,polyCovMATG2R]

    polyfitBTG2L, polyCovMBTG2L =np.polyfit(xVPTG2[0],yTG2[0][1],deg=1,cov=True)
    polyfitBTG2M, polyCovMBTG2M =np.polyfit(xVPTG2[1],yTG2[1][1],deg=1,cov=True)
    polyfitBTG2R, polyCovMBTG2R =np.polyfit(xVPTG2[2],yTG2[2][1],deg=1,cov=True)
    polyfitBTG2 = [polyfitBTG2L,polyfitBTG2M,polyfitBTG2R]
    polyCovMBTG2= [polyCovMBTG2L,polyCovMBTG2M,polyCovMBTG2R]

    polyfitCTG2L, polyCovMCTG2L =np.polyfit(xVPTG2[0],yTG2[0][2],deg=1,cov=True)
    polyfitCTG2M, polyCovMCTG2M =np.polyfit(xVPTG2[1],yTG2[1][2],deg=1,cov=True)
    polyfitCTG2R, polyCovMCTG2R =np.polyfit(xVPTG2[2],yTG2[2][2],deg=1,cov=True)
    polyfitCTG2=[polyfitCTG2L,polyfitCTG2M,polyfitCTG2R]
    polyCovMCTG2 = [polyCovMCTG2L,polyCovMCTG2M,polyCovMCTG2R]
    
    polyfitvectorTG2 = [polyfitATG2,polyfitBTG2,polyfitCTG2]
    polyfitCovMTG2= [polyCovMATG2,polyCovMBTG2,polyCovMCTG2]

    fitvector = [fitvector875, fitvectorTG1, fitvectorTG2]
    linspace = [linspace875, linspaceTG1, linspaceTG2]
    #print '******VERTICAL SCAN*******'
    for pm, posMon in enumerate(fitvector):
        for lm,lossMon in enumerate(posMon):
            #print 'plot['+str(lm)+','+str(pm)+']'
            roots = []
            for sec,section in enumerate(linspace[pm]):
                f7axarr[f][lm,pm].plot(section,lossMon[sec](section),'k-',alpha=0.5)
                if sec>0:
                    root = (lossMon[sec]-lossMon[sec-1]).roots
                    f7axarr[f][lm,pm].plot(root[0],lossMon[sec](root[0]),'o')
                    #print '    root = '+str(root[0])
                    roots.append(root[0])
            center = float(roots[0]+(0.5*(roots[1]-roots[0])))
            centery = float(lossMon[1](center))
            #print '        center = ' + str(center)
            #print type(center)
            #print type(centery)
            annotation_string = r"%0.2f" % (center)
            f7axarr[f][lm,pm].plot(center,centery,'o')
            f7axarr[f][lm,pm].annotate(annotation_string, (center-2, centery-0.005), color='red')

    SlugVFits=[polyfitvector875,polyfitvectorTG1,polyfitvectorTG2]
    SlugVCovM=[polyfitCovM875,polyfitCovMTG1,polyfitCovMTG2]

    f=f+1

#Save fits for analysis:
np.save("TargetFits_FinHFits.npy",FinHFits)
np.save("TargetFits_FinHCovM.npy",FinHCovM)
np.save("TargetFits_SlugHFits.npy",SlugHFits)
np.save("TargetFits_SlugHCovM.npy",SlugHCovM)
np.save("TargetFits_SlugVFits.npy",SlugVFits)
np.save("TargetFits_SlugVCovM.npy",SlugVCovM)
#filename = "TargetFits.pickle"
#with open(filename,'w') as f:
#    pickle.dump(FinHFits,f)
#    pickle.dump(FinHCovM,f)
#    pickle.dump(SlugHFits,f)
#    pickle.dump(SlugHCovM,f)
#    pickle.dump(SlugVFits,f)
#    pickle.dump(SlugVCovM,f)


plt.show()
pp.savefig(fig2)
pp.savefig(fig0[0])
pp.savefig(fig0[1])
pp.savefig(fig0[2])
pp.savefig(fig6[0])
pp.savefig(fig6[1])
pp.savefig(fig7[0])
pp.close()

