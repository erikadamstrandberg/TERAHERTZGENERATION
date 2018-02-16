
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('Agg')
import scipy.constants as const
import matplotlib.pyplot as mplot

def main():
    maxtime = 1000
    maxlength = 500
    dim = [maxtime, maxlength]

    J = np.zeros(dim)
    E = np.ones(dim)
    n = np.ones(dim)
    ny = 1/10

    dt = 1
    dz = 1

    pulselength = 10
    k = 1

    renderstart = 0
    renderstop = 10

    def jtestc(ny, k, t):
        return k/ny  - k/ny*np.exp(-ny*t)
    def jtest(ny, t):
        return np.exp(-ny*t)

    t = np.arange(500)

    Jtestc = np.zeros(maxtime)
    Jtestc[t] = jtestc(ny, k, t)
    
    # for time in t:
    #     Jtestc[time] = jtestc(ny, k, time)

    for t in range(1, J.shape[0]-1):
        for z in range(J.shape[1]):
            J[t, z] = (1 - ny*dt/2)/(1 + ny*dt/2)*J[t-1, z] + dt/2*(n[t-1, z] + n[t, z])*E[t-1, z]/(1 + ny*dt/2)
    mplot.figure(1)
    
    # mplot.subplot(211)
    mplot.plot(J[:, 0], 'b', alpha = 0.4)
    # mplot.plot(Jtest, 'r--')
    mplot.plot(Jtestc, 'g--')
    mplot.ylim(-0.1*10, 10)
    mplot.xlim(0, 100)
    mplot.grid(True)
    mplot.xlabel('$\\frac{t}{\delta t}$')
    mplot.ylabel('$J(t)$')
    
    # mplot.subplot(212)
    # mplot.plot(np.abs(J[:,0]), 'r')
    # mplot.ylim(-1, 1)
    # mplot.xlim(0, J.shape[1])
    # mplot.grid(True)

    filedir = 'jimg/test'
    filename = '.png'
    f = filedir + filename
    mplot.savefig(f)
    mplot.clf()
        

        
if __name__ == '__main__':
    main()
