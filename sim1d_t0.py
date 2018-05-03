import numpy as np
import scipy.constants as const
import matplotlib as matlib
#matlib.use('Agg')                    # This allows the standalone application to plot and save figs.
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
from time import strftime, localtime
from copy import deepcopy
import quickplot

mplot.rcParams['mathtext.fontset'] = 'cm'
mplot.rcParams['font.size'] = 20


def runsim(
        dz = 0.1,               # Space step
        dt = 0.1 - 1e-2,        # Time step
        time = 4000,            # Temporal length of simulation window in units of dt
        size = 4000,           # Spatial length of simulation window in units of dz
        pulsestart = 100,     # Start position of pulse
        pulselength = 13000,    # Spatial length of pulse
        nu = 0,               # Collision rate
        wavelength = 800e-9,  # Laser wavelength
        ramplength = 1/100,    # Electron density initial ramp length as a fraction of size
        rampdamp = 1/10,      # Electron density ramp parameter, (1-exp(-rampdamp*z))
        plasmastopp = deepcopy(size),
        plottime = 0,
        fname = '',
        savesteps = 100,
        t0real = 60e-15):
    

    # Constants
    LIGHTSPEED = const.speed_of_light
    EPSILON = const.epsilon_0

    # Step zero: Establish plasma units. Define omega_0.
    f = LIGHTSPEED/wavelength           # laser frequency
    omegareal = 2*np.pi*f               # laser angular frequency
    OMEGAPRIM = 1                       # this is the normalised omega, use this everywhere in the code
    omega_0 = omegareal/OMEGAPRIM       # this is the laser omega, use this as argument in punits

    #maxlength = 4e-4 # 0.4 mm.
    #maxplasmalength = punit.splasma(maxlength, omega_0)
    #size = int(maxplasmalength/dz)


    


    # Initialise the various fields to be calculated
    size = int(size)
    time = int(time)

    Sample1real = 3e-6
    Sample1 = int(punit.splasma(Sample1real,omega_0))
    Sample2real = 10e-6
    Sample2 = int(punit.splasma(Sample2real,omega_0))
    Sample3real = 100e-6
    Sample3 = int(punit.splasma(Sample3real,omega_0))
    
    pulsesize = punit.sreal(pulselength*dz, omega_0)

    tmin = punit.tplasma((Sample3real + pulsesize)/LIGHTSPEED, omega_0)/dt # ~13000
    tmax = punit.tplasma(punit.sreal(size*dz/2, omega_0)/LIGHTSPEED, omega_0)/dt # ~17500    

    # tmin is the time when the entire pulse has passed s3
    # tmax is the time when the pulse first reaches the edge
    # both expressed in terms of dt
    
    # Time is the arithmetic average of tmax and tmin defined above
    size = int(2/dz * (Sample3 + punit.splasma(pulsesize, omega_0)))
    time = int((Sample3 + size/2 + pulselength)/2*3)
    # This size is slightly longer than the furthest sample point + the size of the pulse.
    # size and time are now defined; we can start initialising the needed matrices.
    dim = [int(time),int(size)]
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

    I0 = 4e18 
    NatREAL = 3e25
    E0REAL = np.sqrt(2*I0/(EPSILON*LIGHTSPEED))
    plasmastart = int(size/2)
    plasmastopp = plasmastart + int(punit.splasma(1e-5, omega_0)/dz)
    E0 = punit.Eplasma(E0REAL,omega_0)
    
    pulselength = int(pulselength/10*4)
    pulsestart = int(size/2) - pulselength
    t0real = 10e-15
    t0 = punit.tplasma(t0real,omega_0)
    Laser.Gauss_forward(E,B,E0,pulselength,pulsestart,OMEGAPRIM,t0,dt,dz) # Sets up the laser pulse in the window
    Natpunit = punit.nplasma(NatREAL,omega_0)
    Nat = np.ones(size)*Natpunit
    Rampfunctions.Ramp_exp(plasmastart,plasmastopp,rampdamp,Natpunit,Nat,size,dt) # Creates a ramp for the electron density
    #mplot.plot(E)
    #mplot.plot(Nat)
    #mplot.savefig('Epresim.png')
    #mplot.axvline(Sample3/dz + plasmastart, linestyle = '--')
    #print(time)
    #mplot.show()
    #return 0
    # These lines plot the initial E-field and the atom density.
    time = 15000
    #mplot.savetxt(Nat)
    #mplot.axis([int(size/2) - pulselength, size, -0.015, 0.015])

    bar = ChargingBar('Simulation running', max = time)
    print(str(datetime.now())+': Beginning simulation.')
    #space = np.arange(size)*dz
    #fig, ax = mplot.subplots()
    #eplot, = ax.plot(space, E, label=r'$E \mathrm{\ [arb.\ u]}$')
    #natplot, = ax.plot(space, Nat, '--', label=r'$n_\mathrm{at}\mathrm{}$')
    #ax.axis([0, len(E)*dz, -2.5e-2, 2.5e-2])
    #mplot.xlabel(r'$\mathrm{Space\ [\delta z]}$')
    #mplot.ylabel(r'$\mathrm{E [arb. u.]}$')
    #ax.legend(loc = 'lower right')
    
    #mplot.axvline(plasmastart + Sample3/dz)
    #fig.savefig('start_2c.pdf', format='pdf', bbox_inches='tight')

    
    #print((Sample3+pulsesize)/dt)
    #print(plasmastart + Sample3/dz)
    Etera1 = np.zeros(time)
    Etera2 = np.zeros(time)
    Etera3 = np.zeros(time)

    simtimes = np.linspace(0.0, 1.0, savesteps)*time
    simtimes = simtimes.astype(int)

    savex = True

    timeinit = strftime('%H%M', localtime())
    
    for i in range(1,time):
        # Calculate all fields for current time
        E = SpaceSolver.E(E,B,J,dt,dz)
        B = SpaceSolver.B(E,B,dt,dz)
        ne = SpaceSolver.N(E,Nat,Ni0,Ni1,Ni2,ne,W1,W2,W3,omega_0,dt)
        J = SpaceSolver.J(E,J,ne,nu,dt,dz)
        E[0] = 0
        B[0] = 0
        
        Etera1[i-1] = E[int(plasmastart+Sample1/dz)]
        Etera2[i-1] = E[int(plasmastart+Sample2/dz)]
        Etera3[i-1] = E[int(plasmastart+Sample3/dz)]
        # This saves the E-field at three different times.
        
        if i in simtimes and savex and False:
            Etot = np.vstack([Etot, E])
            Btot = np.vstack([Btot, B])
            Jtot = np.vstack([Jtot, J])
            Ni0tot.append(Ni0)
            Ni1tot.append(Ni1)
            netot = np.vstack([netot, ne[0]])
        bar.next()
        if i%14000 == 0 :
            mplot.plot(E)
            mplot.plot(netot)
            mplot.axvline(Sample3/dz + plasmastart, linestyle = '--')
            mplot.show()
            mplot.clf()
    netot = ne[0]
    bar.next()
    bar.finish()
    
    print(str(datetime.now())+': Simulation complete.')
    #np.savetxt('ne' + timeinit + '.csv',  netot)

    ne1 = ne[0][int(plasmastart+Sample1/dz)]
    ne2 = ne[0][int(plasmastart+Sample2/dz)]
    ne3 = ne[0][int(plasmastart+Sample3/dz)]
    neE = np.array([ne1,ne2,ne3])
    
    filename = 't0sweep/bm_' + timeinit + 't0' + str(t0real) + '_'
    #np.savetxt(filename + 's' + str(1) + '.csv', Etera1)
    #np.savetxt(filename + 's' + str(2) + '.csv', Etera2)
    np.savetxt(filename + 's' + str(3) + '.csv', Etera3)
    #np.savetxt(filename + 'n.csv', neE)
    #np.savetxt(filename + 'Efinal.csv', E)
    
    #print(netot[0])
    if fname and savex and False:
        #plog('Saving E-field for all time and all space as ' + fname + '.')
        plotnsave(Etot, savetext = True, filename = fname + 'E')
        plotnsave(Btot, savetext = True, filename = fname + 'B')
        plotnsave(Jtot, savetext = True, filename = fname + 'J')
        plotnsave(netot, savetext = True, filename = fname + 'ne')

def plog(msg):
    print(str(datetime.now()) + ': ' + msg)

runsim()
