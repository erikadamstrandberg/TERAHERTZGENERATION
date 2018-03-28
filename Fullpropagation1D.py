import numpy as np
import scipy.constants as const
import matplotlib as matlib
matlib.use('Agg')                    # This allows the standalone application to plot and save figs.
import matplotlib.pyplot as mplot
from pylab import *
import Laser
import Ionization as Ion
import SpaceSolver
import Plasmaunit as punit
import Rampfunctions
import plotnsave
from progress.bar import ChargingBar

from datetime import datetime

from sys import argv

#%% Full propagation! With W1 and W2
def main():
    #print(sys.argv[1][0])
    
    # Define the length and size of the simulation window. Units are in dt and dz respectively.
    TIME = 200
    SIZE = 400
    plott = np.arange(TIME)
    plotz = np.arange(SIZE)
    dim = [TIME,SIZE]
    
    # Initialise the various fields to be calculated
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

    dt = double(1) #Time step
    dz = double(1) #Spatial step. Note: dz = dt is the magic time step in plasma units

    c = const.speed_of_light
    epsilon = const.epsilon_0 
    LAMBDA = 800e-9                     # Wavelength for for the laser pulse
    f = c/LAMBDA
    OMEGAREAL = 2*np.pi*f
    OMEGAPRIM = 1                       # this is the plasma omega, use this everywhere in the code
    OMEGA_0 = OMEGAREAL/OMEGAPRIM       # this is the arbitrary omega, use this as argument in punits
    t0REAL = 50e-16 
    I0 = 4e18 
    NatREAL = 7e26
    E0REAL = np.sqrt(2*I0/(epsilon*c))
    PULSELENGTH = 100                  # Space given to the pulse in the simulation window
    PULSESTART = 100                   # Places the pulse at a certain z
    
    nu = 0                              # Collision freq.
    
    PLASMASTART = PULSELENGTH+PULSESTART    # Were the atom density starts
    PLASMASTOPP = SIZE                      # Were it stops. Note: for PLASMASTOPP = SIZE the plasma extents all the way to the end of the simlation window 
    RAMPLENGTH = 60                        # Length of ramp
    RAMP_DAMP = 0.1                         # 1-exp(-RAMP_DAMP*z) for the exponential ramp

    E0 = punit.Eplasma(E0REAL,OMEGA_0)
    t0 = punit.tplasma(t0REAL,OMEGA_0)
    Laser.Gauss_forward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt,dz) # Sets up the laser pulse in the window
    
    Natpunit = punit.nplasma(NatREAL,OMEGA_0)
    Nat = np.ones(SIZE)*Natpunit
    Rampfunctions.Ramp_exp(PLASMASTART,PLASMASTOPP,RAMP_DAMP,Natpunit,Nat,SIZE,dt) # Sets up the atom density

    #mplot.plot(plotz,Nat)             # Things for checking the setup
    #mplot.plot(plotz,E)
    #Nkritisk = punit.nreal(1,OMEGA_0)
    #mplot.savefig("start")

    #EREALenv = punit.Ereal(E,OMEGA_0)
    #IREAL = epsilon*c*EREALenv**2/2    # This is the real pulse!
    #mplot.plot(plotz,E)                
    #mplot.plot(plotz,Nat)              #Check the setup!
    #Nkritisk = punit.nreal(1,OMEGA_0)  #Check the atom density
    bar = ChargingBar('Simulation running', max = TIME)
    
    print(str(datetime.now())+': Beginning simulation.')
    for i in range(1,TIME):
            
        E = SpaceSolver.E(E,B,J,dt,dz)
        B = SpaceSolver.B(E,B,dt,dz)
        ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,ne,W1,W2,W3,OMEGA_0,dt)
        J = SpaceSolver.J(E,J,ne,nu,dt,dz)

        Etot[i] = E
        Btot[i] = B
        Jtot[i] = J
        Ni0tot[i] = Ni0[0]
        Ni1tot[i] = Ni1[0]
        netot[i] = ne[0]
        bar.next()
        
    print(str(datetime.now())+': Simulation complete.')
    bar.finish()
    z = np.arange(len(Etot[0]))
    plotnsave.plotnsave(z, Etot[140], filename = 'etot_t140')
    mplot.clf()
    
def energy_total_1d(F):
    if np.ndim(F) == 0:
        print('Field seems to be a scalar. Please make sure it\'s a vector.')
        print('Error in energy_total_1d')
        return 0
    elif np.ndim(F) > 1:
        print('Field seems to be a matrix. Please make sure it\'s a vector.')
        print('Error in energy_total_1d')
        return 0
    else:
        return np.sum(F**2)/2
    

if __name__ == '__main__':
    main()
