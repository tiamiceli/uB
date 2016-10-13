#!/usr/local/bin/python

__author__ = 'miceli'

import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import myRootStyle
import scipy.stats as stats
from mpl_toolkits.mplot3d.axes3d import Axes3D
import pandas as pd
import twoDspline

#Force this style on all histograms
mysty = myRootStyle.myRootStyle()
rt.gROOT.SetStyle("mysty")
rt.gROOT.ForceStyle();





#Read in HARP DATA
momentumANDthetaPanda = pd.read_csv('/Users/tiamiceli/Documents/Code/hornOff/HARP_DATA.txt')
momentumANDtheta = momentumANDthetaPanda.as_matrix(columns=momentumANDthetaPanda.columns[1:])


momentumBoundries = np.array([0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 4.00, 5.00, 6.50])
momentumMid = np.zeros(len(momentumBoundries)-1)
for i in range(0,len(momentumBoundries)-1):
    momentumMid[i] = ((momentumBoundries[i] + momentumBoundries[i+1] )/2.)

thetaBoundries = np.array([0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21])
thetaMid = np.zeros(len(thetaBoundries)-1)
for i in range(0,len(thetaBoundries)-1):
    thetaMid[i] = ((thetaBoundries[i] + thetaBoundries[i+1])/2.)

X,Y = np.meshgrid(thetaMid,momentumMid)

XBounds,YBounds = np.meshgrid(thetaBoundries,momentumBoundries)

Z = momentumANDtheta

#error matrix
errorsPanda=pd.read_csv('/Users/tiamiceli/Documents/Code/hornOff/HARP_ERRORS.txt')
errors = errorsPanda.as_matrix(columns=errorsPanda.columns[0:])

errors = errors.reshape(78,78)
errorsT = errors.T

errorsDiag = np.diagonal(errors)
errorsSigma = np.sqrt(errorsDiag)
errorsSigma = errorsSigma.reshape(13,6)

#Make fake data
randseed = 5
#n=2
n = 5000
cv = momentumANDtheta.reshape(1,78)

myfakedataMany, chi2Many = twoDspline.getFakeData(randseed, n, cv, errors)
#print myfakedataMany

myfakedataManylist = []
#
for sim in myfakedataMany:
    myarray = np.array(sim)
    myarray = myarray.reshape(13,6)
    myfakedataManylist.append(myarray)

myfakedataManynps = np.asarray(myfakedataManylist)



#Read in beam MC to be reweighted

momIsPiPlus=[1,1,1]
momMomentum=[1.,2.,3.]
momTheta=[0.05, 0.06, 0.07]
nuE=[0.5,0.4,0.3]
nuEfakeWorlds=[]

for i in range(0,len(momIsPiPlus)):
    cvXsec = twoDspline.xsec(momMomentum[i], momTheta[i], momentumANDtheta, momentumBoundries, thetaBoundries, False)

    thisNuEfakeWorlds = []
    for fakeWorld in myfakedataManynps:
        fakeXsec = twoDspline.xsec(momMomentum[i], momTheta[i], fakeWorld, momentumBoundries, thetaBoundries, False)
        thisNuEfakeWorlds.append(nuE[i]*(cvXsec/fakeXsec))

    nuEfakeWorlds.append(thisNuEfakeWorlds)

nuEfakeWorldsnp=np.array(nuEfakeWorlds)

#plot histogram of nuE for each world
figt = plt.figure(figsize=(8,6))
#plt.set_title(r'$E_{\nu_{\mu}$}')
plt.hist(nuEfakeWorldsnp[:,1],linewidth=0,alpha=1,color='yellow',label=r'world 1')
plt.hist(nuEfakeWorldsnp[:,2],linewidth=0,alpha=0.5,color='darkviolet',label=r'world 2')
plt.legend()



#compute average and standard deviation & fit Gaussians
cv_ave = np.zeros(78)
cv_sig = np.zeros(78)

cv_ave = np.mean(myfakedataManynps, axis=0)
cv_sig = np.std(myfakedataManynps, axis=0)

print "shape of cv_ave = "
print cv_ave.shape

listOfFigSubplots = []
listOfAxesSubplots = []

for p in range(0,len(myfakedataManynps[0])):
    figCombined, axesCombined = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=True)
    #figCombined.set_size_inches(10.5,5.5)
    a = 0
    b = 0
    for t in range(0,len(myfakedataManynps[0][0])):

        title = r'p = ' + str(momentumMid[p]) + r' GeV, $\theta$ = '+ str(thetaMid[t]) + ' Rad'
        axesCombined[a][b].set_title(title, fontsize=10)
        n, bins, patches = axesCombined[a][b].hist(myfakedataManynps[:,p,t],alpha=1,color='mediumturquoise', bins=50)
        #axesCombined[a][b].set_ylim([0,350])
        axesCombined[a][b].tick_params(axis='both', which='major', labelsize=7)
        num = len(myfakedataManynps[:,p,t])*(bins[1]-bins[0])
        gausbins = np.linspace(bins[0],bins[-1],100)
        gaus = stats.norm.pdf(gausbins,loc=cv_ave[p][t],scale=cv_sig[p][t])
        axesCombined[a][b].plot(gausbins,gaus*num)
        b += 1
        if b > 2:
            a += 1
            b = 0

    figCombined.text(0.5, 0.03, r'Differential Cross Section []', ha='center')
    figCombined.text(0.02, 0.5, r'Event Count', va='center', rotation='vertical')
    listOfFigSubplots.append(figCombined)
    listOfAxesSubplots.append(axesCombined)

figchi2 = plt.figure(figsize=(8,6))
n, bins, patches = plt.hist(chi2Many,bins=50,color='orchid',histtype='stepfilled',label=r'multi-universe $\chi^2$')
num = len(chi2Many)*(bins[1]-bins[0])
chi2funcBins=np.linspace(bins[0],bins[-1],200)
chi2func = stats.chi2.pdf(chi2funcBins,df=78)
plt.plot(chi2funcBins,chi2func*num,label=r'$\chi^2$(dof=78)')
plt.legend()
df = stats.chi2.fit(chi2Many,70)
print "fit df = " + str(df)

plt.xlabel(r'$\chi^2$')
plt.ylabel(r'Number of sim worlds per $\chi^2$')
plt.title(r'$\chi^2$')


plt.show()