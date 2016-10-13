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




#make surface
points = 100
myP = []
myT = []
splinedCrossSection = np.zeros((points,points))

for t in range(0,points):
    myT.append(thetaMid[0] + float(t) *(thetaMid[5]-thetaMid[0])/float(points))

for p in range(0,points):
    myP.append(momentumMid[0] + float(p)*(momentumMid[12]-momentumMid[0] )/float(points))

i = 0
for p in myP:

    j = 0
    for t in myT:
        wantPlots = False
        splinedCrossSection[i][j] = twoDspline.xsec(p,t,momentumANDtheta, momentumBoundries,thetaBoundries, wantPlots)
        j+=1
    i+=1

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1,1,1, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, color='lightgreen', alpha=0.5)
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('PiPlus Cross Section (raw)')


myNPT = np.array(myT)
myNPP = np.array(myP)
myNPX,myNPY = np.meshgrid(myNPT,myNPP)
fig2 = plt.figure(figsize=(8,6))
ax2 = fig2.add_subplot(1,1,1, projection='3d')
ax2.plot_surface(myNPX, myNPY, splinedCrossSection, rstride=1, cstride=1, color='lightblue', alpha=0.25)
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('PiPlus Cross Section (splined)')


fig5 = plt.figure(figsize=(8,6))
plt.pcolor(myNPX,myNPY, splinedCrossSection)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('PiPlus Cross Section (splined)')

 #test to see if 1000 worlds create the central value with standard deviation reported by the covariance matrix

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

#compute average and standard deviation
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



print "Results!!!"
for p in range(0,len(momentumMid)):
    for t in range(0,len(thetaMid)):
        print "input xsec["+str(momentumMid[p])+"]["+str(thetaMid[t])+"]: " + str(momentumANDtheta[p][t]) + " +/- " \
            + str(errorsSigma[p][t]) + " ave sim-world xsec = " + str(cv_ave[p][t]) + " +/- " + str(cv_sig[p][t])


plt.show()