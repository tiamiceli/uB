import numpy as np
import ROOT
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from root_numpy import tree2rec
import matplotlib
from scipy import optimize

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


EXT_events  = 62666
NUMI_events = 290771
# POT: [4000,4500] : BNB 1.41E18, NUMI 2.46E18
# POT: [4500,4700] : BNB 3.10E18, NUMI 2.37E18
pot_NUMI = 2.46+2.37

TMIN = 3
TMAX = 23

PEMIN = 50
PEMAX = 800

f_numi = ROOT.TFile("/Users/tiamiceli/Data/TimingNeutrinos/beam_timing_run2_numi.root")
t_numi = f_numi.Get("_tree")
df_NUMI = pd.DataFrame(tree2rec(t_numi,selection='_pe_total > %i & _pe_total < %i & dt<%i & dt>%i'%(PEMIN,PEMAX,TMAX,TMIN)))

f_ext = ROOT.TFile("/Users/tiamiceli/Data/TimingNeutrinos/beam_timing_run2_ext.root")
t_ext = f_ext.Get("_tree")
print f_ext
print t_ext

df_EXT = pd.DataFrame(tree2rec(t_ext,selection='_pe_total > %i & _pe_total < %i & dt<%i & dt>%i'%(PEMIN,PEMAX,TMAX,TMIN)))

dt_NUMI = df_NUMI['_dt'].ravel()
dt_EXT = df_EXT['_dt'].ravel()

# 20 is the usec that we are saving flashed for



# rate in counts / usec
EXT_rate = float(len(dt_EXT)) / ( (TMAX-TMIN) * float(EXT_events) )
# uncertainty on the EXT rate:
EXT_err = np.sqrt(float(len(dt_EXT))) / ( (TMAX-TMIN) * float(EXT_events) )
normed_err = EXT_err/EXT_rate

print 'EXT total number of flashes in %i us window : %i'%(TMAX-TMIN,len(dt_EXT))
print 'EXT rate [counts / us]   = ',EXT_rate
print 'EXT events : ',EXT_events

# average external rate
#avg_EXT = np.average(vals)
#err_EXT = np.sqrt(avg_EXT * EXT_events) / EXT_events
#normed_err = err_EXT/avg_EXT
#print 'the average external rate is ',avg_EXT
#print 'the error on the external rate is ',err_EXT

#tbin_NUMI = 0.5
tbin_NUMI = 0.2
tmin_NUMI = TMIN
tmax_NUMI = TMAX
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
plt.errorbar(bin_centers,vals,xerr=bin_width/2.,yerr=errs,color='r',fmt='o',label='NuMI Run 2 Trigger Data (Beam-On)  [%0.2e events]'%NUMI_events)
plt.axhspan(1-normed_err,1+normed_err,alpha=0.5,color='b',label='Measured Cosmic Rate (Beam-Off)')

#fitting
trapizoidboundsNUMI=([4,   5.5, 0, 10,    15.25],\
                     [5.5, 10,  2, 15.25, 20])
#pNUMI , eNUMI = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, ydata=vals, bounds=trapizoidboundsNUMI)
pNUMI , eNUMI = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, sigma=errs, ydata=vals, p0=[5.5,6.25,1.3,15,15.75], method='lm')
xdNUMI = np.linspace(2.75,23.25,2000)
vertNUMI = [pNUMI[0],pNUMI[1],pNUMI[3],pNUMI[4]]

riseXdiffNUMI  = pNUMI[1]-pNUMI[0]
riseXdiffNUMIe = np.sqrt(eNUMI[1][1]*eNUMI[1][1]+eNUMI[0][0]*eNUMI[0][0] - 2.0*eNUMI[0][1])

fallXdiffNUMI  = pNUMI[4]-pNUMI[3]
fallXdiffNUMIe = np.sqrt(eNUMI[4][4]*eNUMI[4][4]+eNUMI[3][3]*eNUMI[3][3] - 2.0*eNUMI[3][4])

print "riseXdiffNUMI = " + str(riseXdiffNUMI) + " \pm " + str(riseXdiffNUMIe)
print "fallXdiffNUMI = " + str(fallXdiffNUMI) + " \pm " + str(fallXdiffNUMIe)
plt.plot(xdNUMI, trapizoid_bump(xdNUMI, *pNUMI),color='green')
plt.plot(vertNUMI, trapizoid_bump(vertNUMI, *pNUMI),'o',color='green')

plt.xlim([TMIN,TMAX])
plt.legend(numpoints = 1, fontsize=14)
plt.grid()
plt.ylim([0.5,2])
plt.xlabel(r'Time with respect to the NuMI Trigger Time [$\mu$s]')
plt.ylabel(r'Fractional Flash Count per %.02f $\mu$s \newline with respect to Cosmic Background'%tbin_NUMI)
plt.savefig('NUMI_run2.png')
plt.savefig('NUMI_run2.pdf')
plt.show()
