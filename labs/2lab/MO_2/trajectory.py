from matplotlib import colors
from matplotlib.markers import MarkerStyle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri


x, y = [], []
fig, ax = plt.subplots()
plt.grid(True)
val = 0.
try:
     
    file =  open('trajectory.txt', 'r')
    
    coords = [line.split() for line in file]
    for  (c1, c2) in coords:
        x.append(float(c1))
        y.append(float(c2))
    #fig = plt.figure(figsize = (10,10), dpi = 50)
    #ax = fig.add_subplot(1,1,1)
    x1 = np.arange(-0.1, 2)
    y1 = np.arange(-0.1, 2)
    X, Y = np.meshgrid(x1, y1)
    f = (100 * (Y - X) * (Y - X) + (1 - X) * (1 - X))

    ax.contour(X, Y, f, levels=15, linewidths = 1)
    plt.plot(x,y, '-o' )
   
    file =  open('trajectory1.txt', 'r')
    
    coords = [line.split() for line in file]
    for  (c1, c2) in coords:
        x.append(float(c1))
        y.append(float(c2))
    #plt.scatter(x,y, color = 'r')
    triang = []
    k=0
    for i in range(0,int(len(x)/3)):
        triang.append( [k, k+1, k+2])
        k+=3
    triang = tri.Triangulation(x, y, triang)
    #fig1, ax1 = plt.subplots()
    plt.gca().set_aspect('equal', adjustable='box')
    #ax1.triplot(triang, 'bo-', lw=1)
    plt.triplot(triang, 'ro-', lw=0.7, alpha=1,  markersize=1.5)
    plt.savefig(fname='pic_1.png', dpi = 150)

    
except IOError:
    print("Can't open file")
    

ax.set_xlabel("x",fontsize = 14, loc = "right")
ax.set_ylabel("y", fontsize = 14, loc = "top", rotation='horizontal')


plt.show()



