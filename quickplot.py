# -*- coding: utf-8 -*-
from plotnsave import plotnsave
import Plasmaunit as punit
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as mplot

def main():
    c = const.speed_of_light
    k = 2
    name = 'Sample1'
    LAMBDA = 1e-6 
    f = c/LAMBDA
    OMEGAREAL = 2*np.pi*f
    OMEGAPRIM = 1                       # this is the plasma omega, use this everywhere in the code
    OMEGA_0 = OMEGAREAL/OMEGAPRIM
    dt = 0.079

    #Elist = list(np.loadtxt(name + 'E.csv', delimiter = ','))
    #Blist = list(np.loadtxt(name + 'B.csv', delimiter = ','))
    #Jlist = list(np.loadtxt(name + 'J.csv', delimiter = ','))
    #nelist = list(np.loadtxt(name + 'ne.csv', delimiter = ','))
    #print(Elist)
    #print(len(Elist))
    #plotlist = [np.loadtxt]
    #plotlist = nelist[k]
    #print(nelist[k])
    Elist = np.genfromtxt(name + '.csv')
    axis = [0, 80e12 , 1e-6, 1e0]
    Elist = np.abs(np.fft.fft(Elist))
    x = np.fft.fftfreq(len(Elist)) / dt * OMEGA_0
    mplot.plot(x, Elist)
    mplot.xlabel('Hz')
    mplot.ylabel('|E(f)| (plasma units)')
    mplot.yscale('log')
    mplot.axis(axis)
    mplot.savefig(name + '.png')
    print(axis)
    return 0
    #for k in range(len(Elist)):
        #plotlist = [Elist[k], Jlist, nelist]
        #plotlist = [Elist[k], nelist]
        #plotnsave(plotlist, savepic = True, savetext = False, filename = 'quickplot' + str(k) )
        #plotlistfft = [np.fft.fft(Elist[k]), nelist]
        #plotnsave(plotlist, savepic = True, savetext = False, filename = 'quickplotfft' + str(k) )

if __name__ == '__main__':
    main()
