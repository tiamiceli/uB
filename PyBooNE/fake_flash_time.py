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

BNB_events  = 1000000#1000000#1127678 + 795004
EXT_events  = 1000000#36206 + 9551
# POT: [4000,4500] : BNB 1.41E18, NUMI 2.46E18
# POT: [4500,4700] : BNB 3.10E18, NUMI 2.37E18
pot_BNB  = 1.41+3.10

f = ROOT.TFile("/Users/tiamiceli/Data/TimingNeutrinos/fake_data_20ns.root")
fake=True
t = f.Get("opflash_tree_filtered")
df = pd.DataFrame(tree2rec(t))#,selection='pe_total > 50 & dt<23 & dt>0'))
del t

df_BNB  = df.query('trig_word == 2048')
df_EXT  = df.query('trig_word == 512')

dt_BNB  = df_BNB['dt'].ravel()
dt_BNB = dt_BNB[:BNB_events]
print 'dt_BNB length is '+ str(len(dt_BNB))

dt_EXT = df_EXT['dt'].ravel()
dt_EXT = dt_EXT[:EXT_events]
print 'dt_EXT length is '+ str(len(dt_EXT))

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
#tbin_BNB = 0.005
tbin_BNB = 0.01
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
plt.errorbar(bin_centers,vals,xerr=bin_width/2.,yerr=errs,color='r',fmt='o',label='BNB MC  [%0.2e events]'%BNB_events )
plt.axhspan(1-normed_err,1+normed_err,alpha=0.5,color='b',label='Fake Cosmic Rate (Beam-Off)')

#fitting
trapizoidboundsBNB=([2.5, 3.3, 0, 4,   4.85],\
                    [3.3, 4,   2, 4.85, 6  ])

#pBNB , eBNB = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, ydata=vals,p0=[3.2,3.25,1.3,4.6,4.8], method='lm')
pBNB , eBNB = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, ydata=vals, sigma=errs, p0=[3.0,3.3,1.3,4.6,4.9], method='lm')
#pBNB , eBNB = optimize.curve_fit(f=trapizoid_bump, xdata=bin_centers, ydata=vals, bounds=trapizoidboundsBNB)
xdBNB = np.linspace(2,10,1000)
vertBNB = [pBNB[0],pBNB[1],pBNB[3],pBNB[4]]

riseXdiffBNB  = pBNB[1]-pBNB[0]
riseXdiffBNBe = np.sqrt(eBNB[1][1]**2+eBNB[0][0]**2 - 2.0*eBNB[0][1])

fallXdiffBNB  = pBNB[4]-pBNB[3]
fallXdiffBNBe = np.sqrt(eBNB[4][4]*eBNB[4][4]+eBNB[3][3]*eBNB[3][3] - 2.0*eBNB[3][4])

topWidthBNB = pBNB[3]-pBNB[1]
topWidthBNBe = np.sqrt(eBNB[3][3]**2+eBNB[1][1]**2 - 2.0*eBNB[1][3])

bottomWidthBNB = pBNB[4]-pBNB[0]
bottomWidthBNBe = np.sqrt(eBNB[4][4]**2+eBNB[0][0]**2 - 2.0*eBNB[0][4])

p1 = pBNB[0] + 0.5*(pBNB[1]-pBNB[0])
p1e = np.sqrt((eBNB[0][0]) + (0.5**2 * eBNB[1][1]) + (0.5**2 * eBNB[0][0]))

p2 = pBNB[3] + 0.5*(pBNB[4]-pBNB[3])


fwhm = p2 - p1

print "riseXdiffBNB = " + str(riseXdiffBNB) + " \pm " + str(riseXdiffBNBe)
print "fallXdiffBNB = " + str(fallXdiffBNB) + " \pm " + str(fallXdiffBNBe)
print "topWidthBNB = " + str(topWidthBNB) + " \pm " + str(topWidthBNBe)
print "bottomWidthBNB = " + str(bottomWidthBNB) + " \pm " + str(bottomWidthBNBe)
print "fwhm = "+str(fwhm)
print "p1 = " + str(p1) + " \pm " + str(p1e)

print 'Results for ' +str(BNB_events)+' BNB_events and '+str(EXT_events)+' EXT_events:'
print '0 ns jitter'
print 'parameters (x0,x1,y1,x2,x3) = '+ str(pBNB)
print 'error matrix:'
print eBNB

plt.plot(xdBNB, trapizoid_bump(xdBNB, *pBNB),color='green')
plt.plot(vertBNB, trapizoid_bump(vertBNB, *pBNB),'o',color='green')

plt.xlim([2,10])
plt.legend(numpoints = 1, fontsize=14)
plt.grid()
plt.ylim([0.5,2.5])
plt.xlabel(r'Time with respect to the BNB Trigger Time [$\mu$s]')
plt.ylabel(r'Fractional Flash Count per '+str(tbin_BNB)+' $\mu$s \newline with respect to Cosmic Background')
plt.savefig('BNB_fake.png')
plt.savefig('BNB_fake.pdf')
plt.show()



print "done"
