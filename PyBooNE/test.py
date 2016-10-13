import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv('/Users/tiamiceli/Documents/Code/Dirt4/nflux_Dirt3.csv')
#/Users/tiamiceli/Documents/Code/PyBooNE/neutronFlux_last1000.csv
x = np.array(df['xatplane'][(df['xatplane'] < 3000) &(df['xatplane'] > -3000) & (df['yatplane'] < 4000) &(df['yatplane'] > -4000)])
y = np.array(df['yatplane'][(df['xatplane'] < 3000) &(df['xatplane'] > -3000) & (df['yatplane'] < 4000) &(df['yatplane'] > -4000)])
plt.figure()
plt.hist2d(x,y,bins=50)
plt.colorbar()
plt.show()

print "Hello World!"

M = np.array([[1,2,3],[0,4,5],[0,0,6]])
print "M = "
print M

for c in range(0,3):
    for r in range(0,3):
        print M[r][c]

for c in range(0,3):
    for r in range(0,c+1):
        print M[r][c]


