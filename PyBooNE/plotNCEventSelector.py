import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv('/Users/tiamiceli/Documents/Code/MCC5/BNBnu/csvNCEventSelector.csv')

df


x = np.array(df['xatplane'][(df['xatplane'] < 3000) &(df['xatplane'] > -3000) & (df['yatplane'] < 4000) &(df['yatplane'] > -4000)])
y = np.array(df['yatplane'][(df['xatplane'] < 3000) &(df['xatplane'] > -3000) & (df['yatplane'] < 4000) &(df['yatplane'] > -4000)])
plt.figure()
plt.hist2d(x,y,bins=1000)
plt.colorbar()
plt.show()

print "Hello World!"
