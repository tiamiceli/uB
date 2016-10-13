#!/usr/local/bin/python

__author__ = 'miceli'

import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import myRootStyle
import scipy.interpolate as spint
from mpl_toolkits.mplot3d.axes3d import Axes3D
import pandas as pd
import twoDsplineExtended
import twoDspline

#Force this style on all histograms
mysty = myRootStyle.myRootStyle()
rt.gROOT.SetStyle("mysty")
rt.gROOT.ForceStyle()


#Read in HARP DATA
momentumANDthetaPanda = pd.read_csv('/Users/tiamiceli/Documents/Code/hornOff/HARP_DATA.txt')
momentumANDtheta = momentumANDthetaPanda.as_matrix(columns=momentumANDthetaPanda.columns[1:])
print "momentumANDtheta.shape:"
print momentumANDtheta.shape

#Enter the array boundaries
momentumBoundries = np.array([0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 4.00, 5.00, 6.50])
momentumMid = np.zeros(len(momentumBoundries)-1)
for i in range(0,len(momentumBoundries)-1):
    momentumMid[i] = ((momentumBoundries[i] + momentumBoundries[i+1] )/2.)

thetaBoundries = np.array([0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21])
thetaMid = np.zeros(len(thetaBoundries)-1)
for i in range(0,len(thetaBoundries)-1):
    thetaMid[i] = ((thetaBoundries[i] + thetaBoundries[i+1])/2.)

#############################
##   Plot HARP raw data    ##
#############################
fig0 = plt.figure(figsize=(8,6))
plt.pcolor(thetaBoundries,momentumBoundries, momentumANDtheta)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('HARP PiPlus Cross Section (raw)')


#############################
##  Plot HARP spline only  ##
#############################

#Prepare an 2d array to store the cross section information4
thetaBins = np.linspace(thetaMid[0],thetaMid[-1],100)
pBins = np.linspace(momentumMid[0],momentumMid[-1],100)
cvXsec = np.ndarray([len(thetaBins),len(pBins)])

#Do one test to produce splined plots
test = twoDsplineExtended.xsec(momentumBoundries[1], thetaBoundries[1], momentumANDtheta, momentumBoundries, thetaBoundries, True, False)
print "test cross section at momentum = "+str(momentumBoundries[1])+", theta = "+str(thetaBoundries[1])+" is "+str(test)

#Fill cross section 2d array
for ti,t in enumerate(thetaBins):
    for pi,p in enumerate(pBins):
        cvXsec[pi][ti] = twoDsplineExtended.xsec(p, t, momentumANDtheta, momentumBoundries, thetaBoundries, False, False)

#plot histogram of cross section
fig1 = plt.figure(figsize=(2.75,4))
plt.pcolor(thetaBins,pBins, cvXsec)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('PiPlus Cross Section (splined)')

fig2 = plt.figure(figsize=(8,6))
meshAngle,meshMomentum = np.meshgrid(thetaBins,pBins)
ax2 = fig2.gca(projection='3d')
surf2 = ax2.plot_surface(meshAngle,meshMomentum, cvXsec,rstride=5, cstride=5,cmap=cm.rainbow)
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('PiPlus Cross Section (splined)')
#
# #############################
# ##  Plot all phase space   ##
# #############################
#
# #Prepare an 2d array to store the cross section information
# #allThetaBinsBounds = np.linspace(thetaMid[0],thetaMid[-1],50)#
# #allPBinsBounds = np.linspace(momentumMid[0],momentumMid[-1],50)#
# allThetaBinsBounds = np.linspace(0,(np.pi/6.),201)
# allPBinsBounds = np.linspace(0,9,201)
#
# allThetaBinsMid = np.zeros(len(allThetaBinsBounds)-1)
# for i in range(0,len(allThetaBinsBounds)-1):
#     allThetaBinsMid[i] = ((allThetaBinsBounds[i] + allThetaBinsBounds[i+1])/2.)
#
# allPBinsMid = np.zeros(len(allPBinsBounds)-1)
# for i in range(0,len(allPBinsBounds)-1):
#     allPBinsMid[i] = ((allPBinsBounds[i] + allPBinsBounds[i+1])/2.)
#
# cvXsecExtended = np.ndarray([len(allThetaBinsMid),len(allPBinsMid)])
#
# #Fill cross section 2d array
# for ti,t in enumerate(allThetaBinsMid):
#     for pi,p in enumerate(allPBinsMid):
#         cvXsecExtended[pi][ti] = twoDsplineExtended.xsec(p, t, momentumANDtheta, momentumBoundries, thetaBoundries, False, True)
#
# #plot histogram of cross section
# fig3 = plt.figure(figsize=(8,6))
# plt.pcolor(allThetaBinsMid,allPBinsMid, cvXsecExtended)
# plt.colorbar()
# plt.xlabel('Theta (radians)')
# plt.ylabel('Momentum (GeV/c)')
# plt.title('PiPlus Cross Section (splined)')
#
# fig4 = plt.figure(figsize=(8,6))
# meshAngleAll,meshMomentumAll = np.meshgrid(allThetaBinsMid,allPBinsMid)
# ax4 = fig4.gca(projection='3d')
# surf4 = ax4.plot_surface(meshAngleAll,meshMomentumAll, cvXsecExtended,rstride=5, cstride=5,cmap=cm.rainbow)
# plt.xlabel('Theta (radians)')
# plt.ylabel('Momentum (GeV/c)')
# plt.title('PiPlus Cross Section (splined)')
#
#

##############################
## Plot SW Parameterization ##
##############################

#Prepare an 2d array to store the cross section information
#allThetaBinsBounds = np.linspace(thetaMid[0],thetaMid[-1],100)#
#allPBinsBounds = np.linspace(momentumMid[0],momentumMid[-1],100)#
allThetaBinsBounds = np.linspace(0,(np.pi/6.),201)
allPBinsBounds = np.linspace(0,9,201)

allThetaBinsMid = np.zeros(len(allThetaBinsBounds)-1)
for i in range(0,len(allThetaBinsBounds)-1):
    allThetaBinsMid[i] = ((allThetaBinsBounds[i] + allThetaBinsBounds[i+1])/2.)

allPBinsMid = np.zeros(len(allPBinsBounds)-1)
for i in range(0,len(allPBinsBounds)-1):
    allPBinsMid[i] = ((allPBinsBounds[i] + allPBinsBounds[i+1])/2.)

swXsec = np.ndarray([len(allPBinsMid),len(allThetaBinsMid)])

print "swXsec.shape"
print swXsec.shape

#Fill cross section 2d array
for ti,t in enumerate(allThetaBinsMid):
    for pi,p in enumerate(allPBinsMid):
        #swXsec[pi][ti] = twoDsplineExtended.sw(p, t)
        swXsec[pi][ti] = twoDsplineExtended.extsw(p, t)

#plot histogram of cross section
fig5 = plt.figure(figsize=(8,6))
plt.pcolor(allThetaBinsMid,allPBinsMid, swXsec)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
#plt.title('PiPlus SW Cross Section (splined)')
plt.title('PiPlus Extended SW Cross Section (splined)')

fig6 = plt.figure(figsize=(8,6))
meshAngleAll,meshMomentumAll = np.meshgrid(allThetaBinsMid,allPBinsMid)
ax6 = fig6.gca(projection='3d')
surf6 = ax6.plot_surface(meshAngleAll,meshMomentumAll, swXsec,rstride=5, cstride=5,cmap=cm.rainbow)
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('PiPlus SW Cross Section (splined)')

##################################
## Compute volume in eSW & HARP ##
##################################
computedVoleSW = np.ndarray([len(momentumMid),len(thetaMid)])
evaleSW = np.ndarray([len(momentumMid),len(thetaMid)])

volspline = spint.RectBivariateSpline(allPBinsMid,allThetaBinsMid, swXsec)
z=volspline(allPBinsMid,allThetaBinsMid)
for idxlowt,lowt in enumerate(thetaBoundries[:-1]):
    for idxlowp,lowp in enumerate(momentumBoundries[:-1]):

        vol = volspline.integral(lowp,momentumBoundries[idxlowp+1],lowt,thetaBoundries[idxlowt+1])
        xs = vol/( (momentumBoundries[idxlowp+1] - lowp) * (thetaBoundries[idxlowt+1] - lowt) )
        computedVoleSW[idxlowp,idxlowt] = xs
        evaleSW[idxlowp,idxlowt] = volspline.__call__(momentumMid[idxlowp],thetaMid[idxlowt])

computedVoleSW2 = np.ndarray([len(allPBinsMid),len(allThetaBinsMid)])
evaleSW2 = np.ndarray([len(allPBinsMid),len(allThetaBinsMid)])
for idxlowt,lowt in enumerate(allThetaBinsBounds[:-1]):
    for idxlowp,lowp in enumerate(allPBinsBounds[:-1]):
        vol2 = volspline.integral(lowp,allPBinsBounds[idxlowp+1],lowt,allThetaBinsBounds[idxlowt+1])
        xs2 = vol2/( (allPBinsBounds[idxlowp+1] - lowp) * (allThetaBinsBounds[idxlowt+1] - lowt) )
        computedVoleSW2[idxlowp,idxlowt] = xs2
        evaleSW2[idxlowp,idxlowt] = volspline.__call__(allPBinsMid[idxlowp], allThetaBinsMid[idxlowt])

#diff = momentumANDtheta-computedVoleSW
#diff2 = momentumANDtheta-computedVoleSW2
# print "momentumANDtheta"
# print momentumANDtheta
# print "computed vol"
# print computedVoleSW
# #print diff

fig7 = plt.figure(figsize=(8,6))
plt.pcolor(thetaBoundries,momentumBoundries, computedVoleSW)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('Computed spline volumes/area eSW1')

fig75 = plt.figure(figsize=(8,6))
plt.pcolor(thetaBoundries,momentumBoundries, evaleSW)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('Evaluated spline eSW1')

fig76 = plt.figure(figsize=(8,6))
plt.pcolor(thetaBoundries,momentumBoundries, computedVoleSW/evaleSW)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('ratio')

fig8 = plt.figure(figsize=(8,6))
plt.pcolor(allThetaBinsBounds,allPBinsBounds, computedVoleSW2)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('Computed spline volumes/area eSW2')

fig85 = plt.figure(figsize=(8,6))
plt.pcolor(allThetaBinsBounds,allPBinsBounds, computedVoleSW2)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('Evaluated spline eSW2')

#figa = plt.figure(figsize=(8,6))
#plt.volspline

# ##############################
# ## Plot eSW & HARP +/- 1sig ##
# ##############################
# #Read in HARP ERROR
# errorPanda = pd.read_csv('/Users/tiamiceli/Documents/Code/hornOff/HARP_ERRORS.txt')
# errors = errorPanda.as_matrix(columns=errorPanda.columns[1:])
# print "momentumANDtheta.shape:"
# print momentumANDtheta.shape
#
# errors = errors.reshape(78,78)
# errorsT = errors.T
#
# print "errors shape"
# print errors.shape
#
# #np.set_printoptions(threshold='nan')
# #print errorsT.shape
# print "errors!!!!!!="
# for row in range(0,78):
#
#     printrow="[ "
#     for column in range(0,78):
#         temp = str(errors.item((row,column))) + ", "
#         printrow = printrow + temp
#     printrow = printrow + " ],"
#     print printrow
#
# errorsDiag = np.diagonal(errors)
# errorsSigma = np.sqrt(errorsDiag)
# #print errorsSigma
# errorsSigma = errorsSigma.reshape(13,6)
# #print errorsSigma
#
# momentumANDthetaPlusOneSig = momentumANDtheta + errorsSigma
# momentumANDthetaMinusOneSig = momentumANDtheta - errorsSigma
#
# ax.plot_surface(X, Y, momentumANDthetaPlusOneSig, rstride=1, cstride=1, color='yellow', alpha=0.25)
# #cpset = ax.contour(X, Y, momentumANDthetaPlusOneSig, zdir='z', offset=0, cmap='coolwarm')
# #cpset = ax.contour(X, Y, momentumANDthetaPlusOneSig, zdir='x', offset=0.045, cmap='coolwarm')
# #cpset = ax.contour(X, Y, momentumANDthetaPlusOneSig, zdir='y', offset=0.875, cmap='coolwarm')
#
# ax.plot_surface(X, Y, momentumANDthetaMinusOneSig, rstride=1, cstride=1, color='yellow', alpha=0.25)
# #cmset = ax.contour(X, Y, momentumANDthetaMinusOneSig, zdir='z', offset=0, cmap='coolwarm')
# #cmset = ax.contour(X, Y, momentumANDthetaMinusOneSig, zdir='x', offset=0.045, cmap='coolwarm')
# #cmset = ax.contour(X, Y, momentumANDthetaMinusOneSig, zdir='y', offset=0.875, cmap='coolwarm')
#
#

plt.show()