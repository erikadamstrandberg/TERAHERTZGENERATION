import numpy as np
import scipy.constants as const
import matplotlib as matlib
import matplotlib.pyplot as mplot
from pylab import *
import Laser
import Ionization as Ion
import SpaceSolver
import Plasmaunit

#%% Starting parameters! Check your laserpulse! 

TIME = 2000
SIZE = 300
dim = [TIME,SIZE]
E = np.zeros(SIZE)
El = np.zeros(SIZE)
B = np.zeros(SIZE)
J = np.zeros(SIZE)
n = np.zeros(SIZE)
ntemp = np.zeros(SIZE)
Etot = np.zeros(dim)
Btot = np.zeros(dim)
Jtot = np.zeros(dim)
ntot = np.zeros(dim)
k = np.arange(SIZE)
n[100:200] = 1
dt = 1
dz = 1
nu = 0

FREK = 1e20
VAR = 100
Pulselength = 400
(El,OMEGA) = Laser.Gauss(Pulselength,FREK,VAR,dt)
laser = np.arange(len(El))

mplot.plot(laser,El)

#%% Simulate your pulse

for i in range(1,TIME):
    if i < Pulselength:
        E[0] = El[i]
    else:
        E[0] = 0
        
    E = SpaceSolver.E(E,B,J,dt,dz)
    B = SpaceSolver.B(E,B,dt,dz)
    J = SpaceSolver.J(E,J,n,n,nu,dt,dz)
    Etot[i] = E
    Btot[i] = B
    Jtot[i] = J

#%% Blabla

nreal = Plasmaunit.nreal(n[100],OMEGA)
print(nreal)

#%% Plot your pulse

t = 600
mplot.plot(k,Etot[t])
mplot.plot(k,Btot[t])
#mplot.plot(k,J[t])
mplot.plot(k,n)

#%% Test J with n = 0 and E = 0 

SIZE = 10
TIME = 1000
dim = [TIME,SIZE]

J = np.zeros(dim)
E = np.zeros(dim)
n = np.zeros(dim)
ntemp = np.zeros(dim)
v = 1
dt = 0.01
dz = 1
z = np.arange(TIME)

J[0,:] = 1

for i in range(1,len(J[:,0])):
    J[i] = Maxwell.J(E[i-1],J[i-1],n[i-1],ntemp[i-1],v,dt,dz)
    

t = 2
mplot.figure(1)
mplot.plot(z*dt,J[:,t],'r',linewidth=5.0)

z = dt*np.arange(TIME)
J = 1*np.exp(-z*v)
mplot.plot(z,J,'--,b',linewidth=5.0)
mplot.savefig("Jfunkar")

#%% Test J with n and E equal to 1

SIZE = 10
TIME = 1000
dim = [TIME,SIZE]

J = np.zeros(dim)
E = np.ones(dim)
n = np.ones(dim)
ntemp = np.ones(dim)
v = 1
dt = 0.01
dz = 1
z = np.arange(TIME)

J[0,:] = 1.5

for i in range(1,len(J[:,0])):
    J[i] = Maxwell.J(E[i],J[i-1],n[i-1],ntemp[i-1],v,dt,dz)

z = 5
mplot.figure(1)
t = np.arange(TIME)
mplot.plot(t*dt,J[:,z],'r',linewidth=5.0)

z = dt*np.arange(TIME)
J = 1+(1.5-1)*np.exp(-v*z)
mplot.plot(z,J,'--,b',linewidth=5.0)
mplot.savefig("Jfunkar2")

