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
    if not os.path.isfile('neutronFlux_last1000.csv'):
        print 'nflux_last1000.csv not found -- creating'
        allTitle = ['xatplane','yatplane','kErg','event']
        fAll = open('neutronFlux_last1000.csv','a')
        wrA = csv.writer(fAll)
        wrA.writerow(allTitle)
    else:
        fallp = open('neutronFlux_last1000.csv','a')
        wrA = csv.writer(fallp)
        wrA.writerow(line)    

def worker(f):
    try:
        tfile = ROOT.TFile(f,"READ")
        #c = tfile.Get("analysistree/anatree")
        c = tfile.Get("anatree/anatree")
        nEntries = c.GetEntries()
        #if nevts == 0:
        nevts = nEntries
        print 'Found {} events.'.format(nEntries)
        print 'Looping over {} of them.'.format(nevts)
    
        # all proton histograms
        hxplane = []
        hyplane = []
        hKErg   = []
        hEvent  = []
    
        # Loop over events 
        for ientry in xrange(nEntries):
            nb = c.GetEntry(ientry)
            if nb <= 0:
               continue
         
            if ientry >= nevts:
                break
   
            for nentry in xrange(c.geant_list_size):
                if (c.pdg)[nentry] == 2112:
                    if (c.StartPointz[nentry] < -100. and c.EndPointz[nentry] > -100.) or (c.StartPointz[nentry] > -100. and c.EndPointz[nentry] < -100.):
                        zslope = (-100. - c.StartPointz[nentry])/(c.EndPointz[nentry] - c.StartPointz[nentry])
                        hxplane.append(zslope*(c.EndPointx[nentry] - c.StartPointx[nentry]) + c.StartPointx[nentry])
                        hyplane.append(zslope*(c.EndPointy[nentry] - c.StartPointy[nentry]) + c.StartPointy[nentry])
                        hKErg.append(c.Eng[nentry] - c.Mass[nentry])
                        hEvent.append(c.event)
                        
        # write csv
        if len(hxplane) > 0:
            for aline in xrange(len(hxplane)):
                hLine = []
                hLine.append(hxplane[aline])
                hLine.append(hyplane[aline])
                hLine.append(hKErg[aline])
                hLine.append(hEvent[aline])
                kWriter(hLine)
    
    except Exception as e:
        print e
        return False, traceback.format_exc()
    else:
        return True

if __name__ == '__main__':
    sys.exit(main(sys.argv))









