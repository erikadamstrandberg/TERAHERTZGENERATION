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
import matplotlib.ticker as ticker
mplot.rcParams['mathtext.fontset'] = 'cm'


diff = 0.02
START = 0
STOPP = 1+diff
xi = np.arange(START,STOPP,diff)

diff = np.pi/40
START = 0
STOPP = np.pi+diff
phi = np.arange(START,STOPP,diff)

THz_power_xi_phi = np.zeros([len(xi),len(phi)])
Laser_input_xi_phi = np.zeros([len(xi),len(phi)])

for P in range(len(xi)):
    print(P)
    for Q in range(len(phi)):
        print(Q)
        dt = 0.099
        dz = 0.1
        
        c = const.speed_of_light
        epsilon = const.epsilon_0 
        LAMBDA = 800e-9 
        f = c/LAMBDA
        OMEGAREAL = 2*np.pi*f
        OMEGAPRIM = 1                       
        OMEGA_0 = OMEGAREAL/OMEGAPRIM
        t0REAL = 5e-15
        I0 = 1e18
        E0REAL = np.sqrt(2*I0/(epsilon*c))
        NatREAL = 2.7e25
        nuREAL = 100e-15
        nu = punit.omegaplasma(nuREAL,OMEGA_0)
        
        TIME = int(7250/dt)
        SIZE = int(9000/dz)
        PULSELENGTH = int(1250/dz) 
        PULSESTART = int(5500)
        PLASMASTART = PULSELENGTH+PULSESTART
        
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
        Etot = np.zeros(dim)
        netot = np.zeros(dim)
        Z = np.arange(SIZE)*dz
        T = np.arange(TIME)*dt
        Etera = np.zeros(TIME)
        
        threemicroREAL = 3e-6
        threemicro = punit.splasma(threemicroREAL,OMEGA_0)
        threemicro_cut = np.nonzero(Z > threemicro)
        threemicro_i = int(threemicro_cut[0][0])
        
        tenmicroREAL = 10e-6
        tenmicro = punit.splasma(tenmicroREAL,OMEGA_0)
        tenmicro_cut = np.nonzero(Z > tenmicro)
        tenmicro_i = int(tenmicro_cut[0][0])
        
        twentymicroREAL = 20e-6
        twentymicro = punit.splasma(twentymicroREAL,OMEGA_0)
        twentymicro_cut = np.nonzero(Z > twentymicro)
        twentymicro_i = int(twentymicro_cut[0][0])
        
        thirtymicroREAL = 30e-6
        thirtymicro = punit.splasma(thirtymicroREAL,OMEGA_0)
        thirtymicro_cut = np.nonzero(Z > thirtymicro)
        thirtymicro_i = int(thirtymicro_cut[0][0])
        
        PLASMASTOPP = PLASMASTART+threemicro_i
        RAMP_DAMP = 2
        Natpunit = punit.nplasma(NatREAL,OMEGA_0)
        Nat = np.ones(SIZE)*Natpunit
        Rampfunctions.Ramp_exp(PLASMASTART,PLASMASTOPP,RAMP_DAMP,Natpunit,Nat,SIZE,dz)
        
        E0 = punit.Eplasma(E0REAL,OMEGA_0)
        t0 = punit.tplasma(t0REAL,OMEGA_0)
        Laser.TwoC_forward(E,B,E0,xi[P],phi[Q],PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt,dz)
        
        collect_i = PLASMASTART+tenmicro_i
        
        Eft_L = np.fft.fft(E)*2*dt
        Laser_input_xi_phi[P][Q] = sum(np.abs(Eft_L)**2)
        
        for i in range(1,TIME):
            
            E = SpaceSolver.E(E,B,J,dt,dz)
            B = SpaceSolver.B(E,B,dt,dz)
            ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,ne,W1,W2,W3,OMEGA_0,dt)
            J = SpaceSolver.J(E,J,ne,nu,dt,dz)
            E[0] = 0
            B[SIZE-1] = 0
            
            Etera[i-1] = E[collect_i] 
        
        ne_frek = ne[0][collect_i]
        
        N = len(E)
        Efrekreal = np.fft.fftfreq(N)*OMEGA_0/dt
        THzmax = 30e12
        THz_cut = np.nonzero(Efrekreal > THzmax)
        THz_cut_i = int(THz_cut[0][0])
        Eft = np.fft.fft(Etera)*dt*2
        THz_power_xi_phi[P][Q] = sum(np.abs(Eft[0:THz_cut_i])**2)

np.savetxt("Thesis/THz_power_xi_phi_1e18",THz_power_xi_phi,delimiter=",") 
np.savetxt("Thesis/THz_power_xi_phi_1e18_xi",xi,delimiter=",")
np.savetxt("Thesis/THz_power_xi_phi_1e18_phi",phi,delimiter=",")
np.savetxt("Thesis/Laser_input_xi_phi_1e18",Laser_input_xi_phi,delimiter=",")
