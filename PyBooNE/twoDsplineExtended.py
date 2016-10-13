#!/usr/local/bin/python

__author__ = 'miceli'

### This file stores the splining function called xsec and
### the fake data generator getFakeData
###
###

import ROOT as rt
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import myRootStyle
import scipy.stats as stats
from mpl_toolkits.mplot3d.axes3d import Axes3D
import pandas as pd

#Force this style on all histograms
mysty = myRootStyle.myRootStyle()
rt.gROOT.SetStyle("mysty")
rt.gROOT.ForceStyle();



def sw(piMomentum, piAngle, incidentProtonMomentum=8.89):
    a = 220.7
    b = 1.080
    c = 1.000
    d = 1.978
    e = 1.32
    f = 5.572
    g = 0.0868
    h = 9.686
    i = 1.000
    swxsec = a*pow(piMomentum,b)*(1. - ( piMomentum/(incidentProtonMomentum-i) ))*math.exp(-1.*(c*pow(piMomentum,d))/pow(incidentProtonMomentum,e) - f*piAngle*(piMomentum - g*incidentProtonMomentum*pow(np.cos(piAngle),h) ))
    if np.isnan(swxsec):
        swxsec = 0
    elif np.isinf(swxsec):
        swxsec = 1000
    return swxsec

def extsw(piMomentum, piAngle, incidentProtonMomentum=8.89):

    c1 = 5.13;
    c2 = 1.87;
    c3 = 6.67;
    c4 = 1.56;
    c5 = 1.56;
    c6 = 11.9;
    c7 = 0.173;
    c8 = 19.8;
    c9 = 16.0;

    A = c1 - c3*pow(piMomentum,c4)/pow(incidentProtonMomentum,c5) - c6*piAngle*(piMomentum - c7*incidentProtonMomentum*pow(np.cos(piAngle),c8));

    extswxsec = math.exp(A)*pow(piMomentum,c2)*(1 - piMomentum/incidentProtonMomentum)*pow( (1 + piMomentum/incidentProtonMomentum),c9*piAngle*(piMomentum-c7*incidentProtonMomentum*pow(np.cos(piAngle),c8)) );


    if np.isnan(extswxsec):
        extswxsec = 0
    elif extswxsec<0:
        extswxsec = 0
    elif np.isinf(extswxsec):
        extswxsec = 1000

    return extswxsec


def xsec(requestedP, requestedT, PandT, PBoundries, TBoundries, plots=True, extended=False):

    splineConfig=1 #the last bin spacing on the edges of the HARP data are continued into the SW
    #splineConfig=2 #10 bins on all edges of HARP data
    #splineConfig=3 #the finest bin spacing of HARP is used everywhere
    #splineConfig=4 #triple finest bin spacing of HARP for everywhere

    if extended:

        #compute SW outside of HARP data set to put in for "SW extended HARP splines"

        ###############
        ### LOW P   ###
        ###############
        LoPBoundries = np.linspace(0,PBoundries[0],11)
        if splineConfig==2:
            LoPBoundries = np.linspace(0,PBoundries[0],4)
        elif splineConfig==3:
            LoPBoundries = np.linspace(0,PBoundries[0],4)
        elif splineConfig==4:
            LoPBoundries = np.linspace(0,PBoundries[0],10)


        LoPMid = np.zeros(len(LoPBoundries)-1)
        for i in range(0,len(LoPBoundries)-1):
            LoPMid[i] = ((LoPBoundries[i] + LoPBoundries[i+1] )/2.)

        ###############
        ###  HI P   ###
        ###############
        HiPBoundries = np.linspace(PBoundries[-1],8.89,11)
        if splineConfig==2:
            HiPBoundries = np.linspace(PBoundries[-1],9.5,3)
        elif splineConfig==3:
            HiPBoundries = np.linspace(PBoundries[-1],9,11)
        elif splineConfig==4:
            HiPBoundries = np.linspace(PBoundries[-1],9,31)

        HiPMid = np.zeros(len(HiPBoundries)-1)
        for i in range(0,len(HiPBoundries)-1):
            HiPMid[i] = ((HiPBoundries[i] + HiPBoundries[i+1] )/2.)

        ###############
        ### LOW T   ###
        ###############
        LoTBoundries = np.linspace(0,TBoundries[0],11)
        if splineConfig==2:
            LoTBoundries = np.linspace(0,TBoundries[0],2)
        elif splineConfig==3:
            LoTBoundries = np.linspace(0,TBoundries[0],2)
        elif splineConfig==4:
            LoTBoundries = np.linspace(0,TBoundries[0],4)

        LoTMid = np.zeros(len(LoTBoundries)-1)
        for i in range(0,len(LoTBoundries)-1):
            LoTMid[i] = ((LoTBoundries[i] + LoTBoundries[i+1] )/2.)

        ###############
        ###  HI T   ###
        ###############
        HiTBoundries = np.linspace(TBoundries[-1],math.pi,11)
        if splineConfig==2:
            HiTBoundries = np.linspace(TBoundries[-1],0.54,12)
        elif splineConfig==3:
            HiTBoundries = np.linspace(TBoundries[-1],0.54,12)
        elif splineConfig==4:
            HiTBoundries = np.linspace(TBoundries[-1],0.54,34)

        HiTMid = np.zeros(len(HiTBoundries)-1)
        for i in range(0,len(HiTBoundries)-1):
            HiTMid[i] = ((HiTBoundries[i] + HiTBoundries[i+1] )/2.)

        TBoundriesTot = np.append(LoTBoundries[:-1],TBoundries)
        TBoundriesTot = np.append(TBoundriesTot,HiTBoundries[1:])
        TMidTot = np.zeros(len(TBoundriesTot)-1)
        for i in range(0,len(TBoundriesTot)-1):
            TMidTot[i] = ((TBoundriesTot[i] + TBoundriesTot[i+1] )/2.)

        PBoundriesTot = np.append(LoPBoundries[:-1],PBoundries)
        PBoundriesTot = np.append(PBoundriesTot,HiPBoundries[1:])
        PMidTot = np.zeros(len(PBoundriesTot)-1)
        for i in range(0,len(PBoundriesTot)-1):
            PMidTot[i] = ((PBoundriesTot[i] + PBoundriesTot[i+1] )/2.)

        PBoundriesHARP = PBoundries
        PMidHARP = np.zeros(len(PBoundriesHARP)-1)
        for i in range(0,len(PBoundriesHARP)-1):
            PMidHARP[i] = ((PBoundriesHARP[i] + PBoundriesHARP[i+1] )/2.)

        #labelling scheme for other parts of phase space not covered by HARP
        # ^ Momentum
        #__________________
        #|                 |
        #|                 |
        #|     Top         |
        #|                 |
        #|_________________|
        #|L|   HARP   | R  |
        #|_|__________|____|
        #|     Bottom      |
        #|_________________| ->Theta

        # PandTTot = np.ndarray([len(PMidTot),len(TMidTot)])
        # for pi,p in enumerate(PMidTot):
        #     for ti,t in enumerate(TMidTot):
        #         #check for outside HARP
        #         if(p<PBoundries[0] || PBoundries[-1]<p || t<TBoundries[0] || TBoundries[-1]<t):
        #             PandTTot[pi,ti]=sw(p,t)
        #         else:
        #             PandTTot[pi,ti]=




        #compute left sampling
        left = np.ndarray([len(PMidHARP),len(LoTMid)])
        for t,theta in enumerate(LoTMid):
            for p,momentum in enumerate(PMidHARP):
                left[p,t] = extsw(momentum,theta)

        #compute right sampling
        right = np.ndarray([len(PMidHARP),len(HiTMid)])
        for t,theta in enumerate(HiTMid):
            for p,momentum in enumerate(PMidHARP):
                right[p,t] = extsw(momentum,theta)

        #compute top sampling
        top = np.ndarray([len(HiPMid),len(TMidTot)])
        for t,theta in enumerate(TMidTot):
            for p,momentum in enumerate(HiPMid):
                top[p,t] = extsw(momentum,theta)

        #compute bottom sampling
        bottom = np.ndarray([len(LoPMid),len(TMidTot)])
        for t,theta in enumerate(TMidTot):
            for p,momentum in enumerate(LoPMid):
                bottom[p,t] = extsw(momentum,theta)


        #put together the SW outside estimates with the actual HARP data

        LPandT = np.concatenate((left, PandT), axis=1)
        LPandTR = np.concatenate((LPandT,right), axis=1)

        bottomLPandTR = np.concatenate((bottom,  LPandTR), axis=0)
        bottomLPandTRtop = np.concatenate((bottomLPandTR,top),axis=0)

        PandT = bottomLPandTR

        TBoundries = TBoundriesTot
        PBoundries = np.append(LoPBoundries[:-1],PBoundries)
        PBoundries = np.append(PBoundries,HiPBoundries[1:])

        # print "TBoundries:"
        # print TBoundries.shape
        # print TBoundries
        # print "PBoundries:"
        # print PBoundries.shape
        # print PBoundries
        # print "PandT shape:"
        # print PandT.shape
        # print "PandT[0,:]:"
        # print PandT[0,:]
    #select out slices of constant theta and spline them

    splinesConstantT = []
    histConstantT = []

    for t in range(0,len(TBoundries)-1):
        name = "histt_" + str(t)
        bins = len(PBoundries)-1
        titles = "Momentum " + str(TBoundries[t]) + " to " + str(TBoundries[t+1])
        thisHist = rt.TH1F(name,titles,bins, PBoundries)
        thisHist.SetXTitle("Momentum (GeV/c)")
        thisHist.SetYTitle("Entries")
        these_momenta = PandT[:,t]
        #print these_thetas
        p=0
        for element in these_momenta:
            thisbin = p+1
            weightedElement = element #*((TBoundries[t+1] - TBoundries[t]) * (PBoundries[p+1]-PBoundries[p]))
            thisHist.SetBinContent(thisbin,weightedElement)
            p+=1

        histConstantT.append(thisHist)
        thisSpline = rt.TSpline3(thisHist)
        splinename = "splinet_" + str(t)
        splinetitle = "Theta " + str(TBoundries[t]) + " to " + str(TBoundries[t+1]) + ";Momentum (GeV/c);Entries"
        thisSpline.SetNameTitle(splinename,splinetitle)
        splinesConstantT.append(thisSpline)

    #in each theta spline, make a hiso and spline of the value of cross section at p=requestedP
    name = "histp"
    bins = len(TBoundries)-1
    titles = "Momentum at " + str(requestedP)
    histConstantP = rt.TH1F(name,titles,bins,TBoundries)
    histConstantP.SetXTitle("Theta (radians)")
    histConstantP.SetYTitle("Entries")
    t = 0
    for tspline in splinesConstantT:
        thisbin = t+1
        histConstantP.SetBinContent(thisbin,tspline.Eval(requestedP))
        t += 1

    splineConstantP = rt.TSpline3(histConstantP)
    splinename = "splinep"
    splinetitle = "Momementum = " + str(requestedP) + ";Theta (radians);Entries"
    splineConstantP.SetNameTitle(splinename,splinetitle)

    #read the spline over theta at a constant momentum, and read the value for the requested theta:

    crossSection = splineConstantP.Eval(requestedT)

    if plots:
        cConstT = rt.TCanvas( 'cp', 'Splines across momenta', 200, 10, 2000, 2400 )

        cConstT.Divide(3,3)
        for i in range(0,(len(TBoundries)-1)):
            cConstT.cd(i+1)
            splinesConstantT[i].SetLineWidth(2)
            splinesConstantT[i].SetLineColor(rt.kBlue)
            splinesConstantT[i].Draw()
            histConstantT[i].Draw("SAME")

        cConstT.cd(8)
        splineConstantP.SetLineWidth(2)
        splineConstantP.SetLineColor(rt.kRed)
        splineConstantP.Draw()
        histConstantP.Draw("SAME")
        cConstT.Update()
        #cConstT.SaveAs("splines.root")
        cConstT.Print("splines_new.pdf")

    return crossSection

def getFakeData(randseed, n, cv, errmat):
    allFakeData = []
    chi2 = []

    dim = len(errmat[:,1])

    cvArray = rt.TArrayD(dim)
    for c in range(0,dim):
        cvArray[c] = cv[0,c] #CHECK this is weird

    mat = rt.TMatrixD(dim,dim)
    for r in range(0,dim):
        for c in range(0,dim):
            mat[r][c] = errmat[r,c]
    dc = rt.TDecompChol(mat)
    if(not dc.Decompose()):
        print "Can't decompose matrix to generate fake data."
        return allFakeData
    U = dc.GetU()

    rnd = rt.TRandom3()
    rnd.SetSeed(randseed)
    for nWorld in range(0,n):
        fakeData = [0]*dim
        rands = []

        for r in range(0,dim):
            rands.append(rnd.Gaus(0,1))

        for c in range(0,dim):
            rowOfU = 0.
            for r in range(0,c+1):
                rowOfU += U(r,c)*rands[r]
            fakeData[c] = rowOfU + cvArray[c]

            if fakeData[c] < 0 :
                print "Warning: Fake data in bin " + str(c) + " for experiment number " + str(nWorld) \
                      + "is negative. Setting it to 0."
                fakeData[c] = 0.

        allFakeData.append(fakeData)

        fd = np.asarray(fakeData)
        dMinusfd = cv-fd
        invErrMat = np.linalg.inv(errmat)

        thisChi2 = dMinusfd*(invErrMat*dMinusfd.T)
        thisChi2num = float(thisChi2[0][0])


        chi2.append(thisChi2num)

    return allFakeData, chi2