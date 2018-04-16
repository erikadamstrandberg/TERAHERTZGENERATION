# -*- coding: utf-8 -*-
from plotnsave import plotnsave
import Plasmaunit as punit
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as mplot
from os import walk
from plog import plog
from progress.bar import ChargingBar

def plotmatrix(matrix, fname, xmin, xmax, ymin, ymax):
    extent = xmin, xmax, ymin, ymax
    mplot.imshow(matrix, extent = extent, aspect = 'auto')
    mplot.savefig(fname)

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
    
    while True: 
        flaser = omega_0/(2*np.pi)

        input = np.genfromtxt('t0sweep/s3/s3fft.csv')
        input = input[1::, ::]
        maxtime = len(input[0])
        freq_cut = int(3*flaser*maxtime*dt/omega_0)
        maxtimeplot = freq_cut

        fig, ax = mplot.subplots()
        maxes = []
        for i in range(len(input)):
            maxes.append(np.max(input[i]))
        freks = np.fft.fftfreq(maxtime)/dt*omega_0
        '''
        extent = 0, int(maxtimeplot), 1.5e-14, 4e-14
        cax = mplot.imshow(input[::, :int(maxtimeplot):],
                     extent=extent, aspect='auto',
                     vmin=0, vmax=4, cmap='plasma')
        print(np.shape(input[::, :int(maxtimeplot):]))
        print(maxtimeplot)
        fig.colorbar(cax, label=r'|E(\nu)|^2')
        #ax.set_xticks(np.arange(10)*freq_cut/10)
        #ax.set_xlabels(freks[0:maxtimeplot:int(maxtimeplot/10)])
        fig.savefig('t0sweep/s3/s3fft0.png')
        '''
        break
    maxes = np.where(np.array(maxes) < 1)
    # Collect samples
    path = 't0sweep/'
    f = []
    plog('Scanning.')
    for (dirpath, dirnames, filenames) in walk(path):
        f.extend(filenames)
        break
    plog('Sorting.')
    Efin = [file for file in f if '_Efinal' in file]
    s1 = [file for file in f if '_s1' in file]
    s2 = [file for file in f if '_s2' in file]
    s3 = [file for file in f if '_s3' in file]
    maxtime = len(np.genfromtxt(dir + s3[0]))
    plog('Sorted.')
    S3 = np.zeros(maxtime)
    for file in s3:
        data = np.genfromtxt(dir + file)
        data = np.abs(np.fft.fft(data))**2
        S3 = np.vstack([S3, data])

    freq_cut = int(flaser*maxtime*dt/omega_0)
    maxtimeplot = freq_cut*2
    freks = np.arange(10)/10*80e11
    extent = 0, int(maxtimeplot), 1e-15, 10e-15

    fig, ax = mplot.subplots()
    #ax.set_xticks(np.arange(10)) Det här är så man gör
    ax.set_xticks(len(S3[0, ::]))
    #ax.set_xticks(ax.get_xticks/dt/maxtime)
    #ax.set_xticklabels(['cancer'])
    cax = ax.imshow(S3[::, :int(maxtimeplot):], extent=extent,
                    aspect='auto',
                    cmap='plasma')
    fig.colorbar(cax, label=r'|E(\nu)|^2')
    fig.savefig('t0sweep/s3/s3fft02.png')
    plog('Done.')
    return 0

    for file in s3:
        data = np.genfromtxt(dir + file)
        mplot.plot(np.arange(len(data))*dt, data)
        mplot.savefig('t0sweep/test.png')
        mplot.clf()
    for file in Efin:
        data = np.genfromtxt(dir + file)
        mplot.plot(np.arange(len(data))*dt, data)
        mplot.savefig('t0sweep/testfin.png')
    return 0

    for i in range(int(np.max(np.shape(maxes)))):
        print(i)
        Efin[i] = Efin[i][7:-11:]
        print(Efin[int(i)])
        #mplot.(freks, )
    return 0
    #longsim = [file for file in f if 'longsimtestE' in file]
    #Elong = np.loadtxt('t0sweep/longsimtestE.csv', delimiter=',')
    maxtime = len(np.genfromtxt(dir + s1[0]))
    #S1 = unpackdata(s1, dir, maxtime)
    #S2 = unpackdata(s2, dir, maxtime)
    #S3 = unpackdata(s3, dir, maxtime)
    #print(maxtime)

    #np.savetext(dir + 's1/S1.csv', S1)
    #np.savetext(dir + 's2/S2.csv', S2)
    #np.savetext(dir + 's3/S3.csv', S3)
    bar = ChargingBar('Reading data.', max=len(s3))
    np.savetxt(dir + 's3/s3fft.csv', S3)
    # mplot.imshow(plotE, extent = extent, aspect = 'auto')

    return 0
    for thing in s3:
        orig = np.genfromtxt(dir + thing)
        plog('Plotting ' + dir + thing)
        maxtime = len(orig)
        tplasma = np.arange(maxtime)*dt
        treal = punit.treal(tplasma, omega_0)
        #mplot.plot(treal, orig)
        mplot.yscale('log')
        #mplot.savefig(dir + 's3/' + thing[:len(thing)-4] + 'time.png')
        #mplot.clf()
        stuff = np.abs(np.fft.fft(orig))**2
        #stuff = stuff/max(stuff)
        S3 = np.vstack([S3, stuff])
        freks = np.fft.fftfreq(maxtime)/dt*omega_0
        mplot.plot(freks, stuff)
        #mplot.yscale('linear')
        mplot.axis([0, 80e12, 1e-10, max(stuff)*1.2])
        mplot.axvline(1/2/np.pi*omega_0, ls='--')
        mplot.xlabel('Frequency [Hz]')
        mplot.ylabel(r'|E(f)|$^2$ (plasma units)')
        mplot.savefig('s3test_thz.png')
        mplot.axis([0, 3/2/np.pi*omega_0, 1e-10, max(stuff)])
        mplot.savefig('s3test_laser.png')
        break
        if False:
            mplot.title('$t_0$ = 15 fs, z = 100 $\mu$m')
            mplot.savefig(dir + 's3/' + thing[:len(thing)-4] + 'freq_thz.png')
            mplot.axis([0, 2*flaser, 1e-20, max(stuff)])
            mplot.savefig(dir + 's3/' + thing[:len(thing)-4] + 'freq_laser.png')
            #intensity1 = np.vstack([intensity1, stuff])
            mplot.clf()
    plog('Done.')
    return 0
    longspace = np.arange(len(Elong[0]))*dz
    longtimesteps = np.arange(len(Elong))*dt

    z3real = 100e-6
    z3 = int(punit.splasma(z3real, omega_0))
    maxtime = len(np.genfromtxt(dir + s1[0]))
    maxlength = len(np.genfromtxt(dir + Efin[0]))

    longtime = longtimesteps/np.max(longtimesteps)*maxtime
    longtime = longtime.astype(int)

    spacecutmin = 30000
    plotE = Elong[::, spacecutmin::]

    print(np.shape(plotE))

    plotdim = np.shape(plotE)
    plog('Plotting.')

    Elong = np.abs(Elong)**2
    extent = 0, plotdim[0], spacecutmin, plotdim[1] + spacecutmin 
    mplot.imshow(plotE, extent = extent, aspect = 'auto')
    plog('Done!')
    #mplot.contourf(longspace, longtime, Elong)
    mplot.colorbar(label = 'E-field', cmap = 'viridis')
    mplot.ylabel('Space (dz)')
    mplot.xlabel('Time (dt)')
    mplot.savefig('longsim.png')
    return 0
    
    intensity1 = np.zeros([maxtime])
    intensity2 = np.zeros([maxtime])
    intensity3 = np.zeros([maxtime])
    
    tplasma = np.arange(maxtime)*dt
    treal = punit.treal(tplasma, omega_0)
    flaser = (omega_0/2/np.pi)

    
        
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


def unpackdata(input, dir, length):
    output = np.array(length)
    for file in input:
        data = np.genfromtxt(dir + file)
        print(np.shape(data))
        output = np.vstack([output, data])
    return output
    
if __name__ == '__main__':
    main()
