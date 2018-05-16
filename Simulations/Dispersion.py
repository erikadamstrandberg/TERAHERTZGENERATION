import numpy as np
import scipy.constants as const
import matplotlib as matlib
import matplotlib.pyplot as mplot
from pylab import *
import Laser
import SpaceSolver
import Plasmaunit as punit

#%% Full propagation! With W1 and W2

TIME = 8000
SIZE = 15000
plott = np.arange(TIME)
plotz = np.arange(SIZE)
dim = [TIME,SIZE]
E = np.zeros(SIZE)
B = np.zeros(SIZE)
J = np.zeros(SIZE)
W1 = np.zeros(SIZE)
W2 = np.zeros(SIZE)
W3 = np.zeros(SIZE)
Ni2 = (np.zeros(SIZE),np.zeros(SIZE))
Ni1 = (np.zeros(SIZE),np.zeros(SIZE))
Ni0 = (np.ones(SIZE),np.zeros(SIZE))
ne = (np.zeros(SIZE),np.zeros(SIZE))
Ni0tot = np.zeros(dim)
Ni1tot = np.zeros(dim)
netot = np.zeros(dim)
Etot = np.zeros(dim)
Btot = np.zeros(dim)
Jtot = np.zeros(dim)
ntot = np.zeros(dim)
W1tot = np.zeros(dim)

dt = 0.01
dz = 0.01
nu = 0

c = const.speed_of_light
epsilon = const.epsilon_0 
LAMBDA = 800e-9
f = c/LAMBDA
PULSELENGTH = 1000
PULSESTART = 4000
OMEGAREAL = 2*np.pi*f
OMEGAPRIM = 60                      # this is the plasma omega, use this everywhere in the code
OMEGA_0 = OMEGAREAL/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
t0REAL = 30e-15
I0 = 4e18
E0REAL = np.sqrt(2*I0/(epsilon*c))
E0 = punit.Eplasma(E0REAL,OMEGA_0)
t0 = punit.tplasma(t0REAL,OMEGA_0)
Laser.Gauss_forward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt,dz)

mplot.plot(plotz,E)
#mplot.axis([4500,5200,-0.1,0.1])

PLASMASTART = 5000

Omegaplasma = 0.8
ne[0][PLASMASTART:SIZE] = Omegaplasma
ne[1][PLASMASTART:SIZE] = Omegaplasma
mplot.plot(plotz,ne[0])

#%%
for i in range(1,TIME):
    
    E = SpaceSolver.E(E,B,J,dt,dz)
    B = SpaceSolver.B(E,B,dt,dz)
    J = SpaceSolver.J(E,J,ne,nu,dt,dz)
    
    Etot[i] = E
    Btot[i] = B
    Jtot[i] = J
    
#%%
t = 2000
mplot.plot(plotz,Etot[t])
mplot.plot(plotz,ne[0]*1)

#%%

Espace = Etot[400:6000,1000]
Es = np.arange(len(Espace))
mplot.plot(Es,Espace)

Bspace = Btot[400:6000,1000]
Bs = np.arange(len(Bspace))
mplot.plot(Bs,Bspace)

Efs = np.fft.fft(Espace)
Bfs = np.fft.fft(Bspace)
Bfstak = np.conj(Bfs)

Sint = np.multiply(Efs,Bfstak)
Scomplex = (1/2)*sum(Sint)
Sreal = np.real(Scomplex)
print(Sreal)

#%% Calculates fft at Za

t1 = 400
t2 = 3000
Za = 5400
Zb = 5882

N = len(Etot[t1:t2,Za])
Efrek = np.fft.fftfreq(N)
zdiff = (Zb-Za)*dz

T = dt*N
dw = 2*np.pi/T

Espace = Etot[t1:t2,Za]
Es = np.arange(len(Espace))
#mplot.plot(Es,Espace)

Bspace = Btot[t1:t2,Za]
Bs = np.arange(len(Bspace))
#mplot.plot(Bs,Bspace)

Efs = np.fft.fft(Espace)
Bfs = np.fft.fft(Bspace)
mplot.plot(Efrek[0:N/2]*dw*N,np.real(Efs[0:N/2]))
mplot.plot(Efrek[0:N/2]*dw*N,np.imag(Efs[0:N/2]))
mplot.axis([-10,30,-0.001,0.001])
mplot.xlabel("Omega")
mplot.ylabel("E")
mplot.title("Temporal fourier transform at z=za")
mplot.savefig("fft")

Bfstak = np.conj(Bfs)
Sint = np.multiply(Efs,Bfstak)
Scomplex = (1/2)*sum(Sint)
Sreal = np.real(Scomplex)
print(Sreal)

#%% Calculates fft at Zb

Espace2 = Etot[t1:t2,Zb]
Es = np.arange(len(Espace2))
#mplot.plot(Es,Espace2)

Bspace2 = Btot[t1:t2,Zb]
Bs = np.arange(len(Bspace2))
#mplot.plot(Bs,Bspace2)

Efs2 = np.fft.fft(Espace2)
Bfs2 = np.fft.fft(Bspace2)
mplot.plot(Efrek[0:N/2]*dw*N,np.real(Efs2[0:N/2]))
mplot.plot(Efrek[0:N/2]*dw*N,np.imag(Efs2[0:N/2]))
mplot.axis([-10,30,-0.001,0.001])
mplot.xlabel("Omega")
mplot.ylabel("E")
mplot.title("Temporal fourier transform at z=zb")
mplot.savefig("fft2")

Bfstak2 = np.conj(Bfs2)
Sint = np.multiply(Efs2,Bfstak2)
Scomplex = (1/2)*sum(Sint)
Sreal = np.real(Scomplex)
print(Sreal)

#%% The fft extracts the main freq. of the Guassian pulse. OMEGAPRIM

#i = np.argmax(Efsp)
#frek = Efrek[i]*dw*N
#print(frek)

#%%
Efsp = Efs[0:N/2]
Efs2p = Efs2[0:N/2]

Komega = np.divide(Efsp,Efs2p)
phase = np.angle(Komega)
Kunwrap = np.unwrap(phase)/zdiff

Eomegaplot = Efrek[0:N/2]*dw*N
mplot.plot(Kunwrap,Eomegaplot)

k = np.arange(0,300,0.1)
omega = np.sqrt(k**2+Omegaplasma**2)
mplot.plot(k,omega)
x = np.array([0,300])
y = np.array([Omegaplasma,Omegaplasma])
mplot.plot(x,y)
x = np.array([OMEGAPRIM,OMEGAPRIM])
y = np.array([0,200])
mplot.plot(x,y,'black')
mplot.axis([-0.5,3,-0.5,3])
mplot.xlabel("k")
mplot.ylabel("Omega(k)")
mplot.title("Dispersion")
mplot.savefig("Dispersion2")
#np.savetxt("Kunwrap",Kunwrap, delimiter='.')
#np.savetxt("Efrek",Efrekplot, delimiter='.')
