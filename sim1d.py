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
from plotnsave import plotnsave
from progress.bar import ChargingBar
from datetime import datetime
from copy import deepcopy
import quickplot

def runsim(
        dz = 0.1,               # Space step
        dt = 0.1 - 1e-2,               # Time step
        time = 5000,           # Temporal length of simulation window in units of dt
        size = 4000,           # Spatial length of simulation window in units of dz
        pulsestart = 100,     # Start position of pulse
        pulselength = 1000,    # Spatial length of pulse
        nu = 0,               # Collision rate
        wavelength = 800e-9,  # Laser wavelength
        ramplength = 1/100,    # Electron density initial ramp length as a fraction of size
        rampdamp = 1/10,      # Electron density ramp parameter, (1-exp(-rampdamp*z))
        plasmastopp = deepcopy(size),
        plottime = 0,
        fname = '',
        savesteps = 10,
        t0real = 50e-16):
    dim = [int(time),int(size)]

    pulsestart = int(size/2) - pulselength


    
    # Initialise the various fields to be calculated
    size = int(size)
    time = int(time)
    
    #plottime = int(time/2)

    E = np.zeros(size)
    B = np.zeros(size)
    J = np.zeros(size)
    W1 = np.zeros(size)
    W2 = np.zeros(size)
    W3 = np.zeros(size)
    Ni2 = (np.zeros(size),np.zeros(size))
    Ni1 = (np.zeros(size),np.zeros(size))
    Ni0 = (np.ones(size),np.zeros(size))
    ne = np.array([np.zeros(size),np.zeros(size)])
        
    dim = [time,size]
    Ni0tot = []
    Ni1tot = []
    netot = np.zeros(size)
    Etot = np.zeros(size)
    Btot = np.zeros(size)
    Jtot = np.zeros(size)

    LIGHTSPEED = const.speed_of_light
    EPSILON = const.epsilon_0

    f = LIGHTSPEED/wavelength           # laser frequency
    omegareal = 2*np.pi*f               # laser angular frequency
    OMEGAPRIM = 1                       # this is the normalised omega, use this everywhere in the code
    omega_0 = omegareal/OMEGAPRIM       # this is the laser omega, use this as argument in punits
    I0 = 4e18 
    NatREAL = 3e25
    E0REAL = np.sqrt(2*I0/(EPSILON*LIGHTSPEED))
    plasmastart = int(size/2)    # where the electron density starts
    plasmastopp = size                      # Where it stops. Note: for PLASMASTOPP = SIZE the plasma extents all the way to the end of the simlation window
    E0 = punit.Eplasma(E0REAL,omega_0)
    t0 = punit.tplasma(t0real,omega_0)

    Laser.Gauss_forward(E,B,E0,pulselength,pulsestart,OMEGAPRIM,t0,dt,dz) # Sets up the laser pulse in the window
    Natpunit = punit.nplasma(NatREAL,omega_0)
    Nat = np.ones(size)*Natpunit
    Rampfunctions.Ramp_exp(plasmastart,plasmastopp,rampdamp,Natpunit,Nat,size,dt) # Creates a ramp for the electron density
    ne[0] = Nat

    bar = ChargingBar('Simulation running', max = time)
    print(str(datetime.now())+': Beginning simulation.')

    #simtimes = np.linspace(0.0, 1.0, savesteps)*time
    #simtimes = simtimes.astype(int)

    Sample1real = 3e-6
    Sample1 = int(punit.splasma(Sample1real,omega_0))
    
    Sample2real = 10e-6
    Sample2 = int(punit.splasma(Sample2real,omega_0))
    
    Sample3real = 100e-6
    Sample3 = int(punit.splasma(Sample3real,omega_0))
    
    Etera1 = np.zeros(time)
    Etera2 = np.zeros(time)
    Etera3 = np.zeros(time)

    for i in range(1,time):
        # Calculate all fields for current time
        E = SpaceSolver.E(E,B,J,dt,dz)
        B = SpaceSolver.B(E,B,dt,dz)
        ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,ne,W1,W2,W3,omega_0,dt)
        J = SpaceSolver.J(E,J,ne,nu,dt,dz)
        #E[0] = 0
        #B[0] = 0
        #E[len(E)-1] = 0
        #B[len(B)-1] = 0
        #E[0] = 0
        #B[0] = 0

        Etera1[i-1] = E[int(plasmastart+sample1/dz)]
        Etera2[i-1] = E[int(plasmastart+Sample2/dz)]
        Etera3[i-1] = E[int(plasmastart+Sample3/dz)]
        # Save current time
        # if i in simtimes:
        #     Etot = np.vstack([Etot, E])
        #     Btot = np.vstack([Btot, B])
        #     Jtot = np.vstack([Jtot, J])
        #     Ni0tot.append(Ni0)
        #     Ni1tot.append(Ni1)
        #     netot = np.vstack([netot, ne[0]])

        bar.next()
    netot = ne[0]
    bar.next()
    bar.finish()
    
    print(str(datetime.now())+': Simulation complete.')
        
    ne1 = ne[0][int(PLASMASTART+Sample1/dz)]
    ne2 = ne[0][int(PLASMASTART+Sample2/dz)]
    ne3 = ne[0][int(PLASMASTART+Sample3/dz)]
    neE = np.array([ne1,ne2,ne3])

    filename = 'bm_' + timeinit + 't0' + str(t0) + '_'
    plotnsave(Etera1, filename = filename + 's' + str(1))
    plotnsave(Etera2, filename = filename + 's' + str(2))
    plotnsave(Etera3, filename = filename + 's' + str(3))
    plotnsave(neE, filename = filename + 'n')

    
    #z = np.arange(len(Etot[0]))
    plotting = []
    #print(netot[0])
    if fname:
        #plog('Saving E-field for all time and all space as ' + fname + '.')
        plotnsave(Etot, savetext = True, filename = fname + 'E')
        plotnsave(Btot, savetext = True, filename = fname + 'B')
        plotnsave(Jtot, savetext = True, filename = fname + 'J')
        plotnsave(netot, savetext = True, filename = fname + 'ne')

def plog(msg):
    print(str(datetime.now()) + ': ' + msg)
