import numpy as np
import matplotlib as matlib
import matplotlib.pyplot as matplot
import scipy.constants as const
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import Laser
import scipy.fftpack

#%% Single laser pulse

c = const.value("speed of light in vacuum")
SIZE = 500
z = np.arange(SIZE)
E = np.zeros(SIZE)
B = np.zeros(SIZE)
MAXTIME = 500

PULSELENGTH = 400
El = np.zeros(PULSELENGTH)
OMEGA = 1e13
AMPLITUDE = 2
#El = Laser.Laser_gauss_sin(PULSELENGTH,OMEGA,AMPLITUDE)

for i in range(PULSELENGTH):
    El[i] = np.sin(i*OMEGA)

#matplot.plot(z[0:len(El)],El)
for i in range(MAXTIME):
    # BOUNDARY CONDITIONS z=0 in E and z=SIZE in B
    # E[0] = ?
    # B[SIZE-1] = B[SIZE-2]
    for k in range(1,len(E)):
        E[k] = E[k] + (B[k]-B[k-1])
        
    for j in range(len(B)-1):
        B[j] = B[j] + (E[j+1]-E[j])
        
    if i < len(El):
        E[0] = El[i]
    else:
        E[0] = 0

#matplot.plot(z,E,'b')
#matplot.axis([0,SIZE,-1,1])
#matplot.plot(z,B,'g')
Esnok = np.fft.fft(E,SIZE)
Esnoksnygg = np.fft.fftshift(Esnok)
matplot.plot(z,Esnok) 
#%%
#matplot.axis([450,550,-40,80])
Toppar = np.argmax(Esnok)
Efrek = np.fft.fftfreq(SIZE)
Wavefrek = Efrek[Toppar]
print(Wavefrek)


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

#%% Polarized laser!

SIZE = 500
z = np.arange(SIZE)
Ex = np.zeros(SIZE)
Ey = np.zeros(SIZE)
Bx = np.zeros(SIZE)
By = np.zeros(SIZE)
PULSELENGTH = 400
MAXTIME = 500

for i in range(MAXTIME):
    for k in range(len(Ex)):
        Ex[k] = Ex[k] - (By[k]-By[k-1])
        Ey[k] = Ey[k] + (Bx[k]-Bx[k-1])
    for j in range(len(Bx)-1):
        By[j] = By[j] - (Ex[j+1]-Ex[j])
        Bx[j] = Bx[j] + (Ey[j+1]-Ey[j])
    if i < PULSELENGTH:
        xi1 = i/PULSELENGTH
        xi2 = 1-xi1
        Ex[0] = (np.cos(2*np.pi*i*8/200)*np.exp(-(i-PULSELENGTH/2)**2/1e4))*xi1
        Ey[0] = (np.cos(2*np.pi*i*8/200)*np.exp(-(i-PULSELENGTH/2)**2/1e4))*xi2
        
Axes3D.plot(Ex,Ey)
