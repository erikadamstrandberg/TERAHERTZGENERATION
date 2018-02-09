import numpy as np
import scipy.constants as const
import matplotlib as matlib
import matplotlib.pyplot as mplot
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import Laser
import scipy.fftpack
import Ionization as Ion
import j as Jeq

#%% Saves time and space

maxtime = 800
maxlength = 800
dim = [maxtime, maxlength]
E = np.zeros(dim)
B = np.zeros(dim)
J = np.zeros(dim)
n = np.zeros(dim)
ntemp = np.zeros(dim) 
n[:,200:400] = 1
dt = 1
dz = 1
pulselength = 200
z = np.arange(dim[0])
 
for i in range(1,len(E)-1): 
    if i < pulselength:
        E[i,0] = np.cos(2*np.pi*i*8/200)*np.exp(-(i-pulselength/2)**2 /1e4)
    else :
        E[i,0] = 0
    for j in range(1,len(E)):
            E[i, j] = E[i-1, j] - dt/dz * (B[i-1,j] - B[i-1, j-1]) + 0.1*J[i-1,j]
    for j in range(len(B)-1):
            B[i, j] = B[i-1,j] - dt/dz * (E[i, j+1] - E[i, j])
            ntemp[i,j] = n[i,j]
            n[i,j] = n[i,j]
    J[i] = Jeq.Jrad(E[i],J[i],n[i],ntemp[i],2,1,1)
