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
import twoDgaussianProcess
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, Matern, RationalQuadratic

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


# #############################
# ##  Plot HARP smooth only  ##
# #############################

X = thetaMid
y = momentumANDtheta.sum(axis=0)
x = np.linspace(0.03,0.21,100)
X = X.reshape(-1,1)
y = y.reshape(-1,1)
x1 = x.reshape(-1,1)
np.random.seed(1)

##############
## constant kernel + RBF
##############
kernel = ConstantKernel(1.0,(1e-3,1e3)) * RBF(10,(1e-2,1e2))
gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer = 9)
gp.fit(X,y)
y_pred, sigma = gp.predict(x1,return_std=True)
y_pred = y_pred[:,0]
y_pred_minus = y_pred-(1.9600*sigma)
y_pred_plus = y_pred+(1.9600*sigma)

figG = plt.figure()
plt.plot(X, y, 'r.',markersize=10, label="Observations")
plt.plot(x,y_pred, 'k-', label="Prediction")
plt.fill_between(x, y_pred_minus, y_pred_plus, color='green', alpha=0.2, label="95% C.I.")
plt.title('ConstantKernel + RBF')
plt.xlabel('Theta')
plt.ylabel('Cross Section')
plt.legend(loc='upper right')

##############
## constant RBF
##############
kernel = RBF(10,(1e-2,1e2))
gp1 = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer = 9)
gp1.fit(X,y)
y_pred1, sigma1 = gp1.predict(x1,return_std=True)
y_pred1 = y_pred1[:,0]
y_pred_minus1 = y_pred1-(1.9600*sigma1)
y_pred_plus1 = y_pred1+(1.9600*sigma1)

figG1 = plt.figure()
plt.plot(X, y, 'r.',markersize=10, label="Observations")
plt.plot(x,y_pred1, 'k-', label="Prediction")
plt.fill_between(x, y_pred_minus1, y_pred_plus1, color='green', alpha=0.2, label="95% C.I.")
plt.title('RBF')
plt.xlabel('Theta')
plt.ylabel('Cross Section')
plt.legend(loc='upper right')

##############
## rational quadratic
##############
kernel = RationalQuadratic(alpha=1e5, length_scale=0.279)
gp2 = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer = 9)
gp2.fit(X,y)
y_pred2, sigma2 = gp2.predict(x1,return_std=True)
y_pred2 = y_pred2[:,0]
y_pred_minus2 = y_pred2-(1.9600*sigma2)
y_pred_plus2 = y_pred2+(1.9600*sigma2)

figG2 = plt.figure()
plt.plot(X, y, 'r.',markersize=10, label="Observations")
plt.plot(x,y_pred2, 'k-', label="Prediction")
plt.fill_between(x, y_pred_minus2, y_pred_plus2, color='green', alpha=0.2, label="95% C.I.")
plt.title('RationalQuadratic')
plt.xlabel('Theta')
plt.ylabel('Cross Section')
plt.legend(loc='upper right')

##############
## Matern
##############
kernel = Matern(length_scale=0.484, nu=1.5)
gp3 = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer = 9)
gp3.fit(X,y)
y_pred3, sigma3 = gp3.predict(x1,return_std=True)
y_pred3 = y_pred3[:,0]
y_pred_minus3 = y_pred3-(1.9600*sigma3)
y_pred_plus3 = y_pred3+(1.9600*sigma3)

figG3 = plt.figure()
plt.plot(X, y, 'r.',markersize=10, label="Observations")
plt.plot(x,y_pred3, 'k-', label="Prediction")
plt.fill_between(x, y_pred_minus3, y_pred_plus3, color='green', alpha=0.2, label="95% C.I.")
plt.title('Matern Kernel')
plt.xlabel('Theta')
plt.ylabel('Cross Section')
plt.legend(loc='upper right')






#
# #Prepare an 2d array to store the cross section information4
# thetaBins = np.linspace(thetaMid[0],thetaMid[-1],100)
# pBins = np.linspace(momentumMid[0],momentumMid[-1],100)
# #dt = thetaBins[1]-thetaBins[0]
# #dp = pBins[1]-pBins[0]
# cvXsec = np.ndarray([len(thetaBins),len(pBins)])
#
# #Fill cross section 2d array
# for ti,t in enumerate(thetaBins):
#     for pi,p in enumerate(pBins):
#         temp = twoDgaussianProcess.gp2dxsec(p, t, momentumANDtheta, momentumBoundries, thetaBoundries)
#         cvXsec[pi][ti] = temp#*( dt*dp )
#
# #plot histogram of cross section
# fig1 = plt.figure(figsize=(8,6)) #plt.figure(figsize=(2.75,4))
# plt.pcolor(thetaBins,pBins, cvXsec)
# plt.colorbar()
# plt.xlabel('Theta (radians)')
# plt.ylabel('Momentum (GeV/c)')
# plt.title('PiPlus Cross Section (splined)')


# #recreate data bins
# recreatedData = np.ndarray([len(thetaBoundries),len(momentumBoundries)])
# for ti, tmid in enumerate(thetaMid):
#     lo_t = thetaBoundries[ti]
#     hi_t = thetaBoundries[ti+1]
#
#     for pi, pmid in enumerate(momentumMid):
#         lo_p = momentumBoundries[pi]
#         hi_p = momentumBoundries[pi+1]
#
#         recreatedData[pi][ti] =


#
# fig2 = plt.figure(figsize=(8,6))
# meshAngle,meshMomentum = np.meshgrid(thetaBins,pBins)
# ax2 = fig2.gca(projection='3d')
# surf2 = ax2.plot_surface(meshAngle,meshMomentum, cvXsec,rstride=5, cstride=5,cmap=cm.rainbow)
# plt.xlabel('Theta (radians)')
# plt.ylabel('Momentum (GeV/c)')
# plt.title('PiPlus Cross Section (splined)')



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
# print "allThetaBinsMid size = " + str(cvXsecExtended.shape)
#
# #Fill cross section 2d array
# for ti,t in enumerate(allThetaBinsMid):
#     for pi,p in enumerate(allPBinsMid):
#         cvXsecExtended[pi][ti] = twoDsplineExtended.xsec(p, t, momentumANDtheta, momentumBoundries, thetaBoundries, False, True)
#
# #plot histogram of cross section
# fig3 = plt.figure(figsize=(8,6))
# plt.pcolor(allThetaBinsBounds,allPBinsBounds, cvXsecExtended)
# plt.colorbar()
# plt.xlabel('Theta (radians)')
# plt.ylabel('Momentum (GeV/c)')
# plt.title('PiPlus Cross Section Combination Spline')
#
# fig4 = plt.figure(figsize=(8,6))
# meshAngleAll,meshMomentumAll = np.meshgrid(allThetaBinsMid,allPBinsMid)
# ax4 = fig4.gca(projection='3d')
# surf4 = ax4.plot_surface(meshAngleAll,meshMomentumAll, cvXsecExtended,rstride=5, cstride=5,cmap=cm.rainbow)
# plt.xlabel('Theta (radians)')
# plt.ylabel('Momentum (GeV/c)')
# plt.title('PiPlus Cross Section Combination Spline')
#


##############################
## Plot ESW Parameterization ##
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
        swXsec[pi][ti] = twoDgaussianProcess.extsw(p, t)

#plot histogram of cross section
fig5 = plt.figure(figsize=(8,6))
plt.pcolor(allThetaBinsMid,allPBinsMid, swXsec)
plt.colorbar()
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
#plt.title('PiPlus SW Cross Section (splined)')
plt.title('PiPlus ESW Cross Section')

fig6 = plt.figure(figsize=(8,6))
meshAngleAll,meshMomentumAll = np.meshgrid(allThetaBinsMid,allPBinsMid)
ax6 = fig6.gca(projection='3d')
surf6 = ax6.plot_surface(meshAngleAll,meshMomentumAll, swXsec,rstride=5, cstride=5,cmap=cm.rainbow)
plt.xlabel('Theta (radians)')
plt.ylabel('Momentum (GeV/c)')
plt.title('PiPlus ESW Cross Section')


plt.show()