#!/usr/local/bin/python

__author__ = 'miceli'

### This file stores the splining function called xsec and
### the fake data generator getFakeData
###
###

import ROOT as rt
import numpy as np
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


def xsec(requestedP, requestedT, PandT, PBoundries, TBoundries, plots):



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
            thisHist.SetBinContent(thisbin,element)
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

    # if plots:
    #     cConstT = rt.TCanvas( 'cp', 'Splines across momenta', 200, 10, 2000, 2400 )
    #
    #     cConstT.Divide(3,3)
    #     for i in range(0,(len(TBoundries)-1)):
    #         cConstT.cd(i+1)
    #         splinesConstantT[i].SetLineWidth(2)
    #         splinesConstantT[i].SetLineColor(rt.kBlue)
    #         splinesConstantT[i].Draw()
    #         histConstantT[i].Draw("SAME")
    #
    #     cConstT.cd(8)
    #     splineConstantP.SetLineWidth(2)
    #     splineConstantP.SetLineColor(rt.kRed)
    #     splineConstantP.Draw()
    #     histConstantP.Draw("SAME")
    #     cConstT.Update()
    #     #cConstT.SaveAs("splines.root")
    #     cConstT.Print("splines.pdf")

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