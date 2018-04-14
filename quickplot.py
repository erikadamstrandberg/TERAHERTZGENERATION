# -*- coding: utf-8 -*-
from plotnsave import plotnsave
import Plasmaunit as punit
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as mplot
from os import walk
from plog import plog

def main():
    c = const.speed_of_light
    k = 2
    dir = 't0sweep/'
    fname = 'Sample1'
    name = dir + fname
    LAMBDA = 1e-6 
    f = c/LAMBDA
    OMEGAREAL = 2*np.pi*f
    OMEGAPRIM = 1
    omega_0 = OMEGAREAL/OMEGAPRIM
    dz = 0.1
    dt = dz - 1e-2

    # Collect samples
    path = 't0sweep/'
    f = []
    for (dirpath, dirnames, filenames) in walk(path):
        f.extend(filenames)
        break
    Efin = [file for file in f if '_Efinal' in file]
    s1 = [file for file in f if '_s1' in file]
    s2 = [file for file in f if '_s2' in file]
    s3 = [file for file in f if '_s3' in file]
    longsim = [file for file in f if 'longsimtestE' in file]
    Elong = np.loadtxt('t0sweep/longsimtestE.csv', delimiter=',')

    longspace = np.arange(len(Elong[0]))*dz
    longtimesteps = np.arange(len(Elong))*dt

    z3real = 100e-6
    z3 = int(punit.splasma(z3real, omega_0))
    maxtime = len(np.genfromtxt(dir + s1[0]))
    maxlength = len(np.genfromtxt(dir + Efin[0]))

    longtime = longtimesteps/np.max(longtimesteps)*maxtime
    longtime = longtime.astype(int)
    mplot.contourf(longspace, longtime, Elong)
    mplot.colorbar(label = 'E-field')
    mplot.xlabel('Space (dz)')
    mplot.ylabel('Time (dt)')
    mplot.savefig('longsim.png')
    return 0
    
    intensity1 = np.zeros([maxtime])
    intensity2 = np.zeros([maxtime])
    intensity3 = np.zeros([maxtime])

    tplasma = np.arange(maxtime)*dt
    treal = punit.treal(tplasma, omega_0)
    flaser = (omega_0/2/np.pi)

    
    for thing in s3:
        break
        plog('Plotting ' + dir + thing)
        orig = np.genfromtxt(dir + thing)
        maxtime = len(orig)
        tplasma = np.arange(maxtime)*dt
        treal = punit.treal(tplasma, omega_0)
        mplot.plot(treal, orig)
        mplot.yscale('linear')
        mplot.savefig(dir + 's3/' + thing[:len(thing)-4] + 'time.png')
        mplot.clf()
        stuff = np.abs(np.fft.fft(orig))**2
        stuff = stuff/max(stuff)
        freks = np.fft.fftfreq(maxtime)/dt*omega_0
        mplot.plot(freks, stuff)
        mplot.yscale('log')
        mplot.axis([0, 80e12, 1e-20, max(stuff)])
        mplot.axvline(1/2/np.pi * omega_0, ls = '--')
        mplot.xlabel('Frequency [Hz]')
        mplot.ylabel(r'|E(f)|$^2$ (plasma units)')

        mplot.title('$t_0$ = 15 fs, z = 100 $\mu$m')
        mplot.savefig(dir + 's3/' + thing[:len(thing)-4] + 'freq_thz.png')
        mplot.axis([0, 2*flaser, 1e-20, max(stuff)])
        mplot.savefig(dir + 's3/' + thing[:len(thing)-4] + 'freq_laser.png')
        #intensity1 = np.vstack([intensity1, stuff])
        mplot.clf()
        
    for efin in Efin:
        if '1426' in efin:
            plog('Plotting ' + dir + efin)
            stuff = np.loadtxt(dir + efin)
            stuff = np.abs(stuff)
            Nat = np.loadtxt(dir + 'Natlongsim.csv')
            nplot, = mplot.plot(Nat/max(Nat), label = 'Gas density (normalised)')
            maxlength = len(stuff)
            zplasma = np.arange(maxlength)*dz
            zreal = punit.sreal(zplasma, omega_0)
            stuffplot, = mplot.plot(stuff, label = r'$E$-field (arb. units)')
            mplot.yscale('log')
            mplot.legend(handles = [nplot, stuffplot])
            mplot.title('E-field at final time step')
            mplot.ylabel('|E|')
            mplot.xlabel('Length (dz)')
            #mplot.axis([10000, 30000, 1e-11, 1e5])
            
            #mplot.axvline(zreal[int(len(stuff)/2 + z3/dz)])
            mplot.savefig(dir + 'efinal/' + efin[:len(efin)-3] + 'png')
            
            mplot.clf()

    t0values = np.linspace(15e-15, 60e-15, len(s3)+1)
    return 0

    times = np.fft.fftfreq(maxtime)/dt*omega_0
    maxtimes = sum([1 for x in times if x<80e11])

    #print(np.shape(t0values))
    #print(np.shape(times))

    mplot.contourf(times[:len(times)/2 + maxtimes], t0values, intensity1[:, :maxtimes])
    mplot.xlabel('Frequency [Hz]')
    mplot.ylabel('')
    mplot.colorbar()
    mplot.savefig('intensity1.png')
    return 0

    
    if False:
        for thing in Efin:
            stuff = np.genfromtxt(dir + thing)
            mplot.plot(stuff)
            mplot.axvline(int(len(stuff)/2 + z3/dz))
            mplot.savefig(dir + 's3/' +thing[:len(thing)-4] + '.png')
            mplot.clf()
            #break    

    return 0
    #Elist = list(np.loadtxt(name + 'E.csv', delimiter = ','))
    #Blist = list(np.loadtxt(name + 'B.csv', delimiter = ','))
    #Jlist = list(np.loadtxt(name + 'J.csv', delimiter = ','))
    #nelist = list(np.loadtxt(name + 'ne.csv', delimiter = ','))

    Elist = np.genfromtxt(name + '.csv')
    ne = np.genfromtxt('neE.csv')
    axis = [0, 80e12 , 1e-6, 1e-4]
    Elist = np.abs(np.fft.ifft(Elist))/dt
    x = np.fft.fftfreq(len(Elist)) / dt * OMEGA_0
    mplot.plot(x, Elist)
    mplot.xlabel('Hz')
    mplot.ylabel('|E(f)| (plasma units)')
    mplot.yscale('log')
    mplot.axis(axis)
    nyp = 1/(2*np.pi)*np.sqrt(ne[2])*OMEGA_0
    mplot.axvline(nyp, ls = '--')
    mplot.savefig(name + '.png')
    #for k in range(len(Elist)):
        #plotlist = [Elist[k], Jlist, nelist]
        #plotlist = [Elist[k], nelist]
        #plotnsave(plotlist, savepic = True, savetext = False, filename = 'quickplot' + str(k) )
        #plotlistfft = [np.fft.fft(Elist[k]), nelist]
        #plotnsave(plotlist, savepic = True, savetext = False, filename = 'quickplotfft' + str(k) )

if __name__ == '__main__':
    main()
