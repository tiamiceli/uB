import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import cm

gaus = np.random.randn(10000,10000)

plt.hist2d(gaus[0],gaus[1],bins=100,cmap="PuBu_r")

r=[]
g=[]
b=[]
mybins = np.linspace(0,255,num=255)

for i in range(255):
    val = cm.PuBu_r(i)
    r.append(val[0])
    g.append(val[1])
    b.append(val[2])
    #print val
    print str(int(math.floor(r[i]*255))) + " " + str(int(math.floor(g[i]*255))) + " " + str(int(math.floor(b[i]*255)))



plt.figure()
plt.scatter(r,mybins)
plt.scatter(g,mybins)
plt.scatter(b,mybins)

plt.show()