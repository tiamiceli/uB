#!/usr/local/bin/python

__author__ = 'miceli'

import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt
import myRootStyle

neutrinoE = 1. #Average neutrino beam energy GeV
MA = 1.06 #Axial mass GeV/c^2
sinthetaW2 = 0.220 #sin^2(theta_W)

def NCXsec(neutrinoE, MA, q2):
    GF = 0.000011663787 #GeV^-2 Fermi
    MP = 0.938272046 # GeV/c^2proton mass

    SminusU = 4.*MP*neutrinoE - q2

    alpha = 1. - 2.*sinthetaW2
    beta = 1.
    gamma = -(3./2.)*sinthetaW2
    delta = 0.
    tau = q2/(4.*MP*MP)
    kappaP = 1.793
    kappaN = -1.913
    MV = 0.84 # GeV/c^2 is vector dipole mass
    q2overMV2 = q2/(MV*MV)
    gA0 = 1.26 # measured from beta decay

    FV0 = (1.5)*( (kappaP+kappaN)/( (1.+tau)*(1. + q2overMV2)*(1. + q2overMV2) ) )

    FV3 = (0.5)*( (kappaP - kappaN)/( (1. + tau)*(1. + q2overMV2 )*(1. + q2overMV2 ) ))

    GV0 = (1.5)*((1.+kappaN+kappaP)/( (1. + q2overMV2)*(1. + q2overMV2) ))

    GV3 = (0.5)*( (1+kappaP-kappaN)/((1. + q2overMV2)*(1. + q2overMV2)) )

    F2 = alpha*FV3 + gamma*FV0

    F1 = alpha*GV3 + gamma*GV0 -F2

    GA = (0.5)*( gA0/( (1. + (q2/(MA*MA)))*(1. + (q2/(MA*MA))) ) )

    A = (q2/(MP*MP))*( GA*GA*( 1. + tau - F1*F1*( 1. - tau ) + F2*F2*( 1. - tau ) + F1*F2*4.* tau ))

    B = (q2/(MP*MP))*GA*(F1 + F2)

    C = (1./4.)*( GA*GA + F1*F1 +F2*F2*tau)

    crossSection = ((GF*GF*MP*MP)/(np.pi*8*neutrinoE*neutrinoE))\
                   *(A + B*(SminusU/(MP*MP)) + C*(SminusU*SminusU/(MP*MP*MP*MP)) )

    return crossSection

hc = 197.3269718 / 10000000000000000.
hc2 = hc*hc

myq2 = np.array([0.0,0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.3,1.4])

myCrossSections = hc2 * NCXsec(neutrinoE,MA,myq2)

plt.figure()
plt.plot(myq2,myCrossSections)
plt.xlabel(r'$Q^2$ [GeV$^2$/c$^2$]')
plt.ylabel(r'Flux Averaged d$\sigma$/dQ$^2$ [cm$^2$/GeV$^2$/c$^2$]')
plt.title(r'Neutrino-Proton Cross-Section BNL-734')
majorLocatorX = plt.MultipleLocator(0.2)
minorLocatorX = plt.MultipleLocator(0.1)
majorLocatorY = plt.MultipleLocator(0.5e-39)
minorLocatorY = plt.MultipleLocator(0.1e-39)
ax = plt.gca()
ax.xaxis.set_major_locator(majorLocatorX)
ax.xaxis.set_minor_locator(minorLocatorX)
#ax.set_yscale('log')
ax.yaxis.set_major_locator(majorLocatorY)
ax.yaxis.set_minor_locator(minorLocatorY)

plt.show()

myq2paper = np.array([0.45,0.55,0.65,0.75,0.85,0.95,1.05])

myCrossSectionspaper = hc2 * NCXsec(neutrinoE,MA,myq2paper)

print myCrossSectionspaper

