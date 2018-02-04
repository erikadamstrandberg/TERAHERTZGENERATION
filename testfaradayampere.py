import numpy as np
import matplotlib as matlib
import matplotlib.pyplot as matplot
import scipy.constants as const
from pylab import *

#%% Single laser pulse

SIZE = 600
z = np.arange(SIZE)
E = np.zeros(SIZE)
B = np.zeros(SIZE)
PULSELENGTH = 400
MAXTIME = 600

for i in range(MAXTIME):
    for k in range(len(E)):
        E[k] = E[k] + (B[k]-B[k-1])
    for j in range(len(B)-1):
        B[j] = B[j] + (E[j+1]-E[j])
    if i < PULSELENGTH:
        E[0] = np.cos(2*np.pi*i*8/200)*np.exp(-(i-PULSELENGTH/2)**2/1e4)
matplot.plot(z,E)
#matplot.plot(z,B)


#%% Double laser pulse

SIZE = 600
z = np.arange(SIZE)
E = np.zeros(SIZE)
B = np.zeros(SIZE)
PULSELENGTH1 = 400
PULSESTART2 = 100
PULSELENGTH2 = 400
MAXTIME = 600

for i in range(MAXTIME):
    for j in range(len(B)-1):
        B[j] = B[j] + (E[j+1]-E[j])
    for k in range(len(E)):
        E[k] = E[k] + (B[k]-B[k-1])
    if i < PULSELENGTH1:
        E[0] = np.cos(2*np.pi*i*8/200)*np.exp(-(i-PULSELENGTH1/2)**2/1e4)
    elif i > PULSESTART2 and i < PULSELENGTH2+PULSELENGTH2:
        E[0] = np.cos(2*np.pi*i*8/200)*np.exp(-((i-PULSELENGTH2)-PULSELENGTH2/2)**2/1e4)
matplot.plot(z,E)
#matplot.plot(z,B)

#%% Including a J-source term

HALF = 1/2
SIZE = 800
z = np.arange(SIZE)
E = np.zeros(SIZE)
B = np.zeros(SIZE)
J = np.zeros(SIZE)
Ne = np.ones(SIZE)
Ne_temp = 0
Ne[0] = 1
PULSELENGTH = 800
MAXTIME = 400

for i in range(MAXTIME):
    for k in range(len(E)):
        E[k] = E[k] + (B[k]-B[k-1]) - J[k]
    for j in range(len(B)-1):
        B[j] = B[j] + (E[j+1]-E[j])
        Ne_Temp = Ne[j]
        if E[j] > 0.9:
            Ne[j] = j
        else:
            Ne[j] = 0
        J[k] = J[k] + HALF*(Ne_Temp+Ne[j])*E[j]
    if i < PULSELENGTH:
        E[0] = np.cos(2*np.pi*i*8/200)*np.exp(-(i-PULSELENGTH/2)**2/1e4)
matplot.plot(z,E)


#%%

HALF = 1/2
SIZE = 2000
z = np.arange(SIZE)
E = np.zeros(SIZE)
B = np.zeros(SIZE)
J = np.zeros(SIZE)
Ne = np.zeros(SIZE)
Ne[SIZE/2:SIZE-600] = 0.2
PULSELENGTH = 800
MAXTIME = 900

for i in range(MAXTIME):
    for k in range(len(E)):
        E[k] = E[k] + (B[k]-B[k-1]) - J[k]
    for j in range(len(B)-1):
        B[j] = B[j] + (E[j+1]-E[j])
        Ne_Temp = Ne[j]
        J[k] = J[k] + HALF*(Ne_Temp+Ne[j])*E[j]
    if i < PULSELENGTH:
        E[0] = np.cos(2*np.pi*i*8/200)*np.exp(-(i-PULSELENGTH/2)**2/1e4)
matplot.plot(z,E)
matplot.plot(z,Ne)
#matplot.plot(z,B)