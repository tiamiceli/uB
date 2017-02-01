import numpy as np
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

tbin_BNB = 0.15
tmin_BNB = 0.1
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
plt.errorbar(bin_centers,vals,xerr=bin_width/2.,yerr=errs,color='r',fmt='o',label='BNB Trigger Data (Beam-On)  [%.02fE18 POT]'%pot_BNB )
plt.axhspan(1-normed_err,1+normed_err,alpha=0.5,color='b',label='Measured Cosmic Rate (Beam-Off)')
plt.xlim([2,10])
plt.legend(numpoints = 1, fontsize=14)
plt.grid()
plt.ylim([0.9,1.8])
plt.xlabel(r'Time with respect to the BNB Trigger Time [$\mu$s]')
plt.ylabel(r'Fractional Flash Count per 0.15 $\mu$s \newline with respect to Cosmic Background')
plt.savefig('BNB.png')
plt.savefig('BNB.pdf')
plt.show()

tbin_NUMI = 0.5
tmin_NUMI = 0
tmax_NUMI = 25.5
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
plt.errorbar(bin_centers,vals,xerr=bin_width/2.,yerr=errs,color='r',fmt='o',label='NuMI Trigger Data (Beam-On)  [%.02fE18 POT]'%pot_NUMI)
plt.axhspan(1-normed_err,1+normed_err,alpha=0.5,color='b',label='Measured Cosmic Rate (Beam-Off)')
plt.xlim([2.75,23.25])
plt.legend(numpoints = 1, fontsize=14)
plt.grid()
plt.ylim([0.9,1.5])
plt.xlabel(r'Time with respect to the NuMI Trigger Time [$\mu$s]')
plt.ylabel(r'Fractional Flash Count per 0.5 $\mu$s \newline with respect to Cosmic Background')
plt.savefig('NUMI.png')
plt.savefig('NUMI.pdf')
plt.show()
