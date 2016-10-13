#!/usr/local/bin/python


import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

wave = pd.read_csv('/Users/tiamiceli/Documents/Code/RWM/RWMdata.txt')

ts1 = wave.loc[1,'w4000':'w4003']
print ts1

word0 = int(ts1[0])
word1 = int(ts1[1])
word2 = int(ts1[2])
word3 = int(ts1[3])

hex0 = hex(word0)
hex1 = hex(word1)
hex2 = hex(word2)
hex3 = hex(word3)

print "\nhex words"
print hex0
print hex1
print hex2
print hex3

#1 high 16 bits of seconds since 1/1/1970
#2 low 16 bits of seconds since 1/1/1970
#3 high 16 bits of microsec within that second
#4 low 16 bits of microsec within that second

seconds = (word0<<16)|(word1)
usecond = (word2<<16)|(word3)

print "\nseconds:"
print seconds
print "hex seconds:"
print hex(seconds)
print "usecond:"
print usecond
print "hex usecond"
print hex(usecond)
print " "


print time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(seconds))

#unsigned long long integerPieces[4];
#       for(int j=0; j<4; ++j){
#       integerPieces[j] = ACNET_RWM_ARRAY[NSAMPLES+j];
#       }
#
#       unsigned long long my64bitTimeStamp;
#       my64bitTimeStamp =  (integerPieces[0] << 48)
#       |(integerPieces[0] << 32)
#       |(integerPieces[0] << 16)
#       |(integerPieces[0]      );

#print list(wave.columns.values)

#print wave.tail(3)

#numSamples = np.linspace(0,4004,num=4005)
#for num in numSamples:
#    print "w"+str(int(num))+",",