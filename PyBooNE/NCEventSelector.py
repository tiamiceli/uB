#!/usr/local/bin/python

import os
import sys
import inspect
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle
import sys
import csv
import multiprocessing
import traceback

class PyTree(ROOT.TTree):
    def Hi(self):
        return 'Hi'

def main(argv):

    inputFile = ''
    inputFileList = ''
    nevts = 0
    # Parse arguments
    args = argv[1:]
    while len(args) > 0:
        if args[0] == '-s' or args[0] == '--source' :
            if len(args) > 1:
                inputFile = args[1]
                del args[0:2]
        elif args[0] == '-S' or args[0] == '--source-list' :
            if len(args) > 1:
                inputFileList = args[1]
                del args[0:2]
        elif args[0] == '-n' or args[0] == '--nevts' :
            if len(args) > 1:
                nevts = int(args[1])
                del args[0:2]

    if inputFile == '' and inputFileList == '':
        print 'No input file(s) specified. Use "-s" or "--source" to specify one. Additionally, a file list can be supplied with "-S" or "--source-list".'
        return 0

    plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
   
    # Read in file(s) and create the TChain c
    inputFiles = []
    #c = ROOT.TChain("analysistree/anatree") 
    if inputFile != '':
        inputFiles.append(inputFile)
    elif inputFileList != '':
        inputFiles = open(inputFileList).read().splitlines()
    
    pool = multiprocessing.Pool(processes=6)
    for f in inputFiles:
        print 'Adding input file: {}'.format(f)
        pool.apply_async(worker,args=(f,))
    pool.close()
    pool.join()

def kWriter(line):
    if not os.path.isfile('csvNCEventSelector.csv'):
        print 'csvNCEventSelector.csv not found -- creating'
        allTitle = ['subrun','event','pdg','KE','Px','Py','Pz','StartPointx','StartPointy','StartPointz','EndPointx','EndPointy','EndPointz']
        fAll = open('csvNCEventSelector.csv','a')
        wrA = csv.writer(fAll)
        wrA.writerow(allTitle)
    else:
        fallp = open('csvNCEventSelector.csv','a')
        wrA = csv.writer(fallp)
        wrA.writerow(line)    

def worker(f):
    try:
        tfile = ROOT.TFile(f,"READ")
        c = tfile.Get("analysistree/anatree")
        #c = tfile.Get("anatree/anatree")
        #c = ROOT.TChain("anatree/anatree")
        #print type(c)
        #c.Add(f)
        nEntries = c.GetEntries()
        #if nevts == 0:
        nevts = nEntries
        print 'Found {} events.'.format(nEntries)
        print 'Looping over {} of them.'.format(nevts)
    
        # all histograms
        hsubrun = []
        hevent = []
        hpdg = []
        hKE = []
        hPx = []
        hPy = []
        hPz = []
        hStartPointx = []
        hStartPointy = []
        hStartPointz = []
        hEndPointx = []
        hEndPointy = []
        hEndPointz = []

        # Loop over events 
        for ientry in xrange(nEntries):
            nb = c.GetEntry(ientry)
            if nb <= 0:
               continue
         
            if ientry >= nevts:
                break

            if c.ccnc_truth != 1: # if not NC
                continue
            if c.mode_truth != 0: # if not Elastic or Quasi-Elastic
                continue
            if c.hitnuc_truth != 2112 and c.hitnuc_truth != 2212: # if not hit neutron or proton
                continue

            for nentry in xrange(c.geant_list_size):
                mass = 0
                if c.pdg[nentry] == 2112: #if neutron
                    mass = 0.9395654
                if c.pdg[nentry] == 2212: #if proton
                    mass = 0.9382720
                if (c.pdg[nentry] == 2112) or (c.pdg[nentry] == 2212): #if neutron or proton
                    #check that start point is in detector fiducial region
                    if ( 0 < c.StartPointx[nentry]  and c.EndPointx[nentry] < 256.35) \
                            and ( -116.5 < c.StartPointy[nentry] and c.EndPointy[nentry] < 116.5) \
                            and ( 0 < c.StartPointz[nentry] and c.StartPointz[nentry] < 1036.8):

                        if c.process_primary[nentry] != 1: #if particle is not primary
                            continue

                        hsubrun.append(c.subrun)
                        hevent.append(c.event)
                        hpdg.append(c.pdg[nentry])
                        hKE.append(c.Eng[nentry] - mass)
                        hPx.append(c.Px[nentry])
                        hPy.append(c.Py[nentry])
                        hPz.append(c.Pz[nentry])
                        hStartPointx.append(c.StartPointx[nentry])
                        hStartPointy.append(c.StartPointy[nentry])
                        hStartPointz.append(c.StartPointz[nentry])
                        hEndPointx.append(c.EndPointx[nentry])
                        hEndPointy.append(c.EndPointy[nentry])
                        hEndPointz.append(c.EndPointz[nentry])

        # write csv
        if len(hPx) > 0:
            for aline in xrange(len(hPx)):
                hLine = []
                hLine.append(hsubrun[aline])
                hLine.append(hevent[aline])
                hLine.append(hpdg[aline])
                hLine.append(hKE[aline])
                hLine.append(hPx[aline])
                hLine.append(hPy[aline])
                hLine.append(hPz[aline])
                hLine.append(hStartPointx[aline])
                hLine.append(hStartPointy[aline])
                hLine.append(hStartPointz[aline])
                hLine.append(hEndPointx[aline])
                hLine.append(hEndPointy[aline])
                hLine.append(hEndPointz[aline])
                kWriter(hLine)
    
    except Exception as e:
        print e
        return False, traceback.format_exc()
    else:
        return True

if __name__ == '__main__':
    sys.exit(main(sys.argv))









