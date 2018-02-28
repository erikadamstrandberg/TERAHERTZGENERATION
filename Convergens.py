import numpy as np
import scipy.constants as const
import matplotlib as matlib
import matplotlib.pyplot as mplot
from pylab import *
import Laser
import SpaceSolver
import Plasmaunit as punit
import Rampfunctions

#%% Full propagation! With W1 and W2

#n = np.arange(0,1,0.01)
n = [0,0.1]
matdim = [len(n),3]
matmitten = np.zeros(matdim)
matvid75 = np.zeros(matdim)
matvid875 = np.zeros(matdim)

for x in range(len(n)):
    dt = 1-(n[x]*0.9)
    dz = 1-(n[x]*0.9)
    nu = 0
    
    TIME = int(2000+(n[x]*20000))
    SIZE = int(TIME)
    PULSELENGTH = int(SIZE/2)
    RAMPLENGTH = int(SIZE/10)
    
    plott = np.arange(TIME)
    plotz = np.arange(SIZE)
    dim = [TIME,SIZE]
    E = np.zeros(SIZE)
    Etemp = np.zeros(SIZE)
    B = np.zeros(SIZE)
    J = np.zeros(SIZE)
    W1 = np.zeros(SIZE)
    W2 = np.zeros(SIZE)
    W3 = np.zeros(SIZE)
    Ni2 = np.zeros(SIZE)
    Ni2temp = np.zeros(SIZE)
    Ni1 = np.zeros(SIZE)
    Ni1temp = np.zeros(SIZE)
    Ni0 = np.ones(SIZE)
    Ni0temp = np.zeros(SIZE)
    ne = np.zeros(SIZE)
    netemp = np.zeros(SIZE)
    Ni0tot = np.zeros(dim)
    Ni1tot = np.zeros(dim)
    netot = np.zeros(dim)
    Etot = np.zeros(dim)
    Btot = np.zeros(dim)
    Jtot = np.zeros(dim)
    ntot = np.zeros(dim)
    W1tot = np.zeros(dim)
    
    c = const.speed_of_light
    epsilon = const.epsilon_0 
    LAMBDA = 800e-9 
    f = c/LAMBDA
    PULSESTART = 0
    OMEGAREAL = 2*np.pi*f
    OMEGAPRIM = 1                       # this is the plasma omega, use this everywhere in the code
    OMEGA_0 = OMEGAREAL/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
    t0REAL = 50e-15 
    I0 = 4e18 
    NatREAL = 3e25
    E0REAL = np.sqrt(2*I0/(epsilon*c))
    
    E0 = punit.Eplasma(E0REAL,OMEGA_0)
    t0 = punit.tplasma(t0REAL,OMEGA_0)
    Laser.Gauss_forward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt)
    
    #EREALenv = punit.Ereal(E,OMEGA_0)
    #IREAL = epsilon*c*EREALenv**2/2    # This is the real pulse!
    
    #mplot.plot(plotz,E)
    
    PLASMASTART = PULSELENGTH+PULSESTART
    PLASMASTOPP = SIZE
    
    Natpunit = punit.nplasma(NatREAL,OMEGA_0)
    Nat = np.ones(SIZE)*Natpunit
    
    Rampfunctions.Ramp_linear(RAMPLENGTH,PLASMASTART,PLASMASTOPP,Natpunit,Nat,SIZE,dt)
    
    #mplot.plot(plotz,Nat)
    Nkritisk = punit.nreal(1,OMEGA_0)

    for i in range(1,TIME):
        E = SpaceSolver.E(E,B,J,dt,dz)
        B = SpaceSolver.B(E,B,dt,dz)   
        ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,Ni0temp,Ni1temp,ne,W1,W2,W3,OMEGA_0,dt)
        Ni0temp = Ni0
        Ni1temp = Ni1
        J = SpaceSolver.J(E,J,ne,netemp,nu,dt,dz)
        netemp = ne
    
    matmitten[x,0] = E[int(SIZE/2)]
    matmitten[x,1] = J[int(SIZE/2)]
    matmitten[x,2] = ne[int(SIZE/2)]
    matvid75[x,0] = E[int(SIZE/2+SIZE/4)]
    matvid75[x,1] = J[int(SIZE/2+SIZE/4)]
    matvid75[x,2] = ne[int(SIZE/2+SIZE/4)]
    matvid875[x,0] = E[int(SIZE/2+SIZE/4+SIZE/8)]
    matvid875[x,1] = J[int(SIZE/2+SIZE/4+SIZE/8)]
    matvid875[x,2] = ne[int(SIZE/2+SIZE/4+SIZE/8)]
    
np.savetxt('matmitten.txt', matmitten, delimiter='.')
np.savetxt('matvid75.txt', matvid75, delimiter='.')
np.savetxt('matvid875.txt', matvid875, delimiter='.')
