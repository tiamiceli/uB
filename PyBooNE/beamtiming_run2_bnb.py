import numpy as np
import ROOT
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from root_numpy import tree2rec
import matplotlib
import sys
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

BNB_events  = 1127678 + 795004
EXT_events  = 36206 + 9551
NUMI_events = 276205 + 91004
# POT: [4000,4500] : BNB 1.41E18, NUMI 2.46E18
# POT: [4500,4700] : BNB 3.10E18, NUMI 2.37E18


#f = ROOT.TFile("beam_timing.root")
f = ROOT.TFile("/Users/tiamiceli/Data/TimingNeutrinos/beam_timing_run2_bnb.root")
t = f.Get("_tree")
df = pd.DataFrame(tree2rec(t))#,selection='_pe_total > 50 and _dt<23 and _dt>0'))

print df.shape

good_runs = [8317,8318,8319,8320,8321,8322,8324,8325,8329,8330,8331,8333,8336,8338,8346,8347,8348,8349,8350,8351,8352,8354,8355,8356,8357,8358,8383,8385,8386,8388,8393,8394]

good_string = '(_run == %i) '%good_runs[0]
for n in xrange(1,len(good_runs)):
    good_string += ' or (_run == %i)'%(good_runs[n])

print good_string

#df = df.query(good_string)

#df = df.query('(_run == 8317) or (_run == 8318)  or (_run == 8319) or (_run == 8320) or (_run == 8321) or (_run == 8322) or (_run == 8324) or (_run == 8325) or (_run == 8329) or (_run == 8330) or (_run == 8331) or (_run == 8333) or (_run == 8336) or (_run == 8338)')

print df.shape

df_BNB  = df.query('_trig_word == 2048')
df_NUMI = df.query('_trig_word == 4096')
df_EXT  = df.query('_trig_word == 512')

print 'BNB shape : ',df_BNB.shape
print 'EXT shape : ',df_EXT.shape



df = df.query('(_dt > 0) and (_dt < 23) and (_pe_total > 50) and (_pe_total < 2000)')

print 'DF shape : ',df.shape
#print df
#del t

df_BNB  = df.query('_trig_word == 2048')
df_NUMI = df.query('_trig_word == 4096')
df_EXT  = df.query('_trig_word == 512')

BNB_events = len(df_BNB.groupby(['_event','_run']))
EXT_events = len(df_EXT.groupby(['_event','_run']))
NUMI_events = len(df_NUMI.groupby(['_event','_run']))

#print 'BNB shape : ',df_BNB.shape
#print 'EXT shape : ',df_EXT.shape

#BNB_events = 321568#len(df_BNB.groupby(['_event','_run']))
#EXT_events = 5421#len(df_EXT.groupby(['_event','_run']))
#NUMI_events = len(df_NUMI.groupby(['_event','_run']))

#print 'BNB events : %i'%BNB_events
#print 'EXT events : %i'%EXT_events


pot_BNB  = BNB_events * 0.000003
pot_NUMI = 2.46+2.37



dt_BNB  = df_BNB['_dt'].ravel()
dt_NUMI = df_NUMI['_dt'].ravel()
dt_EXT = df_EXT['_dt'].ravel()

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
tmin_BNB = 2
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
plt.errorbar(bin_centers,vals,xerr=bin_width/2.,yerr=errs,color='r',fmt='o',label='BNB Run 2 Trigger Data (Beam-On)  [%0.2e events]'%BNB_events )
plt.axhspan(1-normed_err,1+normed_err,alpha=0.5,color='b',label='Measured Cosmic Rate (Beam-Off)')

#fitting
#trapizoidboundsBNB=([2.5, 3.225, 0, 4,   4.75],\
#                    [3.225, 4,   2, 4.75, 6  ])

trapizoidboundsBNB=([2.5, 2.5, 0, 4,   4],\
                    [4, 4,   2, 6, 6  ])


pBNB , eBNB = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, ydata=vals,p0=[3.2,3.25,1.3,4.6,4.8], method='lm')#, bounds=trapizoidboundsBNB)
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

plt.xlim([tmin_BNB,tmax_BNB])
plt.legend(numpoints = 1, fontsize=14)
plt.grid()
plt.ylim([0.5,2.5])
plt.xlabel(r'Time with respect to the BNB Trigger Time [$\mu$s]')
plt.ylabel(r'Fractional Flash Count per %.02f $\mu$s \newline with respect to Cosmic Background'%(tbin_BNB))
plt.savefig('BNB_run2.png')
plt.savefig('BNB_run2.pdf')
plt.show()

#sys.ext(0)
