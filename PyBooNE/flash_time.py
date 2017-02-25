import numpy as np
import scipy as sp
from scipy import optimize
import ROOT
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from root_numpy import tree2rec
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# [4000,4500] : 1127678 BNB events, 276205 NuMI events, 36206 cosmic (EXT) events
# 795004 BNB events, 91004 NuMI events, 9551 cosmic (EXT) events

def trapizoid_bump(x, x0, x1, y1, x2, x3):
    return np.piecewise(x, [x < x0,\
                            (x >= x0) & (x < x1),\
                            (x >= x1) & (x < x2),\
                            (x >= x2) & (x < x3),\
                            x>= x3
                            ],\
                        [lambda x:1,\
                         lambda x:x*((y1-1)/(x1-x0)) + (1-x0*((y1 -1)/(x1-x0))),\
                         lambda x:y1,\
                         lambda x:x*((1-y1)/(x3-x2)) + (1-x3*((1-y1)/(x3-x2))),\
                         lambda x:1\
                         ])

BNB_events  = 1127678 + 795004
EXT_events  = 36206 + 9551
NUMI_events = 276205 + 91004
# POT: [4000,4500] : BNB 1.41E18, NUMI 2.46E18
# POT: [4500,4700] : BNB 3.10E18, NUMI 2.37E18
pot_BNB  = 1.41+3.10
pot_NUMI = 2.46+2.37

f = ROOT.TFile("/Users/tiamiceli/Data/TimingNeutrinos/daq_hist_small.root")
t = f.Get("opflash_tree_filtered")
df = pd.DataFrame(tree2rec(t,selection='pe_total > 50 & dt<23 & dt>0'))
del t

df_BNB  = df.query('trig_word == 2048')
df_NUMI = df.query('trig_word == 4096')
df_EXT  = df.query('trig_word == 512')

dt_BNB  = df_BNB['dt'].ravel()
dt_NUMI = df_NUMI['dt'].ravel()
dt_EXT = df_EXT['dt'].ravel()

# 20 is the usec that we are saving flashed for

# rate in counts / usec
EXT_rate = float(len(dt_EXT)) / ( 23. * float(EXT_events) )
# uncertainty on the EXT rate:
EXT_err = np.sqrt(float(len(dt_EXT))) / ( 23. * float(EXT_events) )
normed_err = EXT_err/EXT_rate

print 'EXT total number of flashes in 23 us window : ',len(dt_EXT)
print 'EXT rate [counts / us]   = ',EXT_rate
print 'EXT events : ',EXT_events

# average external rate
#avg_EXT = np.average(vals)
#err_EXT = np.sqrt(avg_EXT * EXT_events) / EXT_events
#normed_err = err_EXT/avg_EXT
#print 'the average external rate is ',avg_EXT
#print 'the error on the external rate is ',err_EXT

#tbin_BNB = 0.15
tbin_BNB = 0.02
tmin_BNB = 2#0.1
tmax_BNB = 10
bins_BNB = np.linspace(tmin_BNB,tmax_BNB, int((tmax_BNB-tmin_BNB)/tbin_BNB) )
vals, bins = np.histogram(dt_BNB,bins=bins_BNB)
vals = vals.astype(float)
bin_width = bins[1]-bins[0]
bin_centers = 0.5*(bins[1:]+bins[:-1])
errs = np.sqrt(vals)
vals /= ( tbin_BNB * float(BNB_events) )
print 'BNB rate [counts / us ]  = ',vals[-3]
vals /= EXT_rate
errs /= ( tbin_BNB * float(BNB_events) )
errs /= EXT_rate
fig = plt.figure(figsize=(10,6))
plt.errorbar(bin_centers,vals,xerr=bin_width/2.,yerr=errs,color='r',fmt='o',label='BNB Run 1 Trigger Data (Beam-On)  [%0.2e events]'%BNB_events )
plt.axhspan(1-normed_err,1+normed_err,alpha=0.5,color='b',label='Measured Cosmic Rate (Beam-Off)')

#fitting
trapizoidboundsBNB=([2.5, 3.3, 0, 4,   4.85],\
                    [3.3, 4,   2, 4.85, 6  ])

pBNB , eBNB = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, ydata=vals,p0=[3.2,3.25,1.3,4.6,4.8], method='lm')
#pBNB , eBNB = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, ydata=vals, bounds=trapizoidboundsBNB)
xdBNB = np.linspace(2,10,1000)
vertBNB = [pBNB[0],pBNB[1],pBNB[3],pBNB[4]]

print "eBNB="
print eBNB

riseXdiffBNB  = pBNB[1]-pBNB[0]
riseXdiffBNBe = np.sqrt(eBNB[1][1]**2+eBNB[0][0]**2 - 2.0*eBNB[0][1])

fallXdiffBNB  = pBNB[4]-pBNB[3]
fallXdiffBNBe = np.sqrt(eBNB[4][4]*eBNB[4][4]+eBNB[3][3]*eBNB[3][3] - 2.0*eBNB[3][4])

print "riseXdiffBNB = " + str(riseXdiffBNB) + " \pm " + str(riseXdiffBNBe)
print "fallXdiffBNB = " + str(fallXdiffBNB) + " \pm " + str(fallXdiffBNBe)

plt.plot(xdBNB, trapizoid_bump(xdBNB, *pBNB),color='green')
plt.plot(vertBNB, trapizoid_bump(vertBNB, *pBNB),'o',color='green')

plt.xlim([2,10])
plt.legend(numpoints = 1, fontsize=14)
plt.grid()
plt.ylim([0.5,2.5])
plt.xlabel(r'Time with respect to the BNB Trigger Time [$\mu$s]')
plt.ylabel(r'Fractional Flash Count per '+str(tbin_BNB)+' $\mu$s \newline with respect to Cosmic Background')
plt.savefig('BNB_run1.png')
plt.savefig('BNB_run1.pdf')
plt.show()

#NUMI

#tbin_NUMI = 0.5
#tbin_NUMI = 0.08
tbin_NUMI = 0.1
tmin_NUMI = 3#0
tmax_NUMI = 23 #25.5
bins_NUMI = np.linspace(tmin_NUMI,tmax_NUMI, int((tmax_NUMI-tmin_NUMI)/tbin_NUMI) )
vals, bins = np.histogram(dt_NUMI,bins=bins_NUMI)
bin_width = bins[1]-bins[0]
vals = vals.astype(float)
bin_centers = 0.5*(bins[1:]+bins[:-1])
errs = np.sqrt(vals)
vals /= ( tbin_NUMI * float(NUMI_events) )
print 'NuMI rate [counts / us ] = ',vals[-3]
vals /= EXT_rate
errs /= ( tbin_NUMI * float(NUMI_events) )
errs /= EXT_rate
fig = plt.figure(figsize=(10,6))
plt.errorbar(bin_centers,vals,xerr=bin_width/2.,yerr=errs,color='r',fmt='o',label='NuMI Run 1 Trigger Data (Beam-On)  [%0.2e events]'%NUMI_events)
plt.axhspan(1-normed_err,1+normed_err,alpha=0.5,color='b',label='Measured Cosmic Rate (Beam-Off)')

#fitting
trapizoidboundsNUMI=([4,   5.8, 0, 10,    15.3],\
                      [5.8, 10,  2, 15.3, 20])
#pNUMI , eNUMI = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, ydata=vals, bounds=trapizoidboundsNUMI)
pNUMI , eNUMI = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, ydata=vals, p0=[5.5,6.25,1.3,15,15.75])
xdNUMI = np.linspace(2.75,23.25,1000)
vertNUMI = [pNUMI[0],pNUMI[1],pNUMI[3],pNUMI[4]]

riseXdiffNUMI  = pNUMI[1]-pNUMI[0]
riseXdiffNUMIe = np.sqrt(eNUMI[1][1]*eNUMI[1][1]+eNUMI[0][0]*eNUMI[0][0] - 2.0*eNUMI[0][1])

fallXdiffNUMI  = pNUMI[4]-pNUMI[3]
fallXdiffNUMIe = np.sqrt(eNUMI[4][4]*eNUMI[4][4]+eNUMI[3][3]*eNUMI[3][3] - 2.0*eNUMI[3][4])

print "riseXdiffNUMI = " + str(riseXdiffNUMI) + " \pm " + str(riseXdiffNUMIe)
print "fallXdiffNUMI = " + str(fallXdiffNUMI) + " \pm " + str(fallXdiffNUMIe)
plt.plot(xdNUMI, trapizoid_bump(xdNUMI, *pNUMI),color='green')
plt.plot(vertNUMI, trapizoid_bump(vertNUMI, *pNUMI),'o',color='green')

#plt.xlim([2.75,23.25])
plt.xlim([3,23])
plt.legend(numpoints = 1, fontsize=14)
plt.grid()
plt.ylim([0.5,2])
plt.xlabel(r'Time with respect to the NuMI Trigger Time [$\mu$s]')
plt.ylabel(r'Fractional Flash Count per '+str(tbin_NUMI)+' $\mu$s \newline with respect to Cosmic Background')
plt.savefig('NUMI_run1.png')
plt.savefig('NUMI_run1.pdf')
plt.show()

print "done"
