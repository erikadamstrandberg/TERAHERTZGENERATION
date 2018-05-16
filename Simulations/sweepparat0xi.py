#%%
import numpy as np
import scipy.constants as const
import matplotlib as matlib
import matplotlib.pyplot as mplot
from pylab import *
import Laser
import SpaceSolver
import Plasmaunit as punit
import Rampfunctions
import Ionization as Ion

dt = 0.0495
dz = 0.05
nu = 0

TIME = int(1250/dt)
SIZE = int(2500/dz)
PULSELENGTH = int(1250/dz)
PULSESTART = int(500/dz)
PLASMASTART = PULSELENGTH + PULSESTART

diff = 0.02
START = 0
STOPP = 1+diff
xi = np.arange(START,STOPP,diff)

diff = 1e-15
START = 0.5e-15
STOPP = 20e-15+diff
t0REAL = np.arange(START,STOPP,diff)

THz_power_xi_t0 = np.zeros([len(xi),len(t0REAL)])
Laser_input_xi_t0 = np.zeros([len(xi),len(t0REAL)])
Ratio_xi_t0 = np.zeros([len(xi),len(t0REAL)])

for P in range(len(xi)):
    for Q in range(len(t0REAL)):
        dim = [TIME,SIZE]
        plott = np.arange(TIME)
        plotz = np.arange(SIZE)
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
        Z = np.arange(SIZE)*dz
        T = np.arange(TIME)*dt
        Etera = np.zeros(TIME)
        
        c = const.speed_of_light
        epsilon = const.epsilon_0 
        LAMBDA = 800e-9 
        f = c/LAMBDA
        OMEGAREAL = 2*np.pi*f
        OMEGAPRIM = 1                       # this is the plasma omega, use this everywhere in the code
        OMEGA_0 = OMEGAREAL/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
        NatREAL = 3e25
        I0 = 4e18
        E0REAL = np.sqrt(2*I0/(epsilon*c))
        
        E0 = punit.Eplasma(E0REAL,OMEGA_0)
        t0 = punit.tplasma(t0REAL[Q],OMEGA_0)
        
        phi = np.pi/2
        Laser.TwoC_forward(E,B,E0,xi[P],phi,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt,dz)
            
        Natpunit = punit.nplasma(NatREAL,OMEGA_0)
        Nat = np.ones(SIZE)*Natpunit
        
        PLASMASTOPP = SIZE
        RAMP_DAMP = 2
        
        Rampfunctions.Ramp_exp(PLASMASTART,PLASMASTOPP,RAMP_DAMP,Natpunit,Nat,SIZE,dz)
        Nkritisk = punit.nreal(1,OMEGA_0)
        
        Sample1real = 3e-6
        Sample1 = int(punit.splasma(Sample1real,OMEGA_0))
        
        Eft_L = np.fft.fft(E)*2*dt
        Laser_input_xi_t0[P][Q] = sum(np.abs(Eft_L)**2)
    
        for i in range(1,TIME):
    
            E = SpaceSolver.E(E,B,J,dt,dz)
            B = SpaceSolver.B(E,B,dt,dz)
            ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,ne,W1,W2,W3,OMEGA_0,dt)
            J = SpaceSolver.J(E,J,ne,nu,dt,dz)
            E[0] = 0
            B[SIZE-1] = 0
            
            Etera[i-1] = E[int(PLASMASTART+Sample1/dz)]
             
        N = len(Etera)
        Efrekreal = np.fft.fftfreq(N)*OMEGA_0/dt
        Eft = np.fft.fft(Etera)*2*dt
        
        THzmax = 30e12
        THz_cut = np.nonzero(Efrekreal > THzmax)
        THz_cut_i = int(THz_cut[0][0])
        THz_power_xi_t0[P][Q] = sum(np.abs(Eft[0:THz_cut_i])**2)
        Ratio_xi_t0[P][Q] = THz_power_xi_t0[P][Q]/Laser_input_xi_t0[P][Q]
        
np.savetxt("THz_power_xit0",THz_power_xi_t0,delimiter=",") 
np.savetxt("Ratio_xi_t0",Ratio_xi_t0,delimiter=",")
np.savetxt("Laser_input_xi_t0",Laser_input_xi_t0,delimiter=",")
np.savetxt("THz_power_xit0_xi",xi,delimiter=",")
np.savetxt("THz_power_xit0_t0",t0REAL,delimiter=",")