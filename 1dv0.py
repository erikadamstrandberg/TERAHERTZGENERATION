# -*- coding: utf-8 -*-
#

import numpy as np
import matplotlib
matplotlib.use('Agg')
import scipy.constants as const
import matplotlib.pyplot as pyplot

def main():
    maxtime = 1000
    maxlength = 500
    dim = [maxtime, maxlength]
    E = np.zeros(dim)
    B = np.zeros(dim)
    J = np.zeros(dim)
    n = np.ones(dim)
    dt = 1
    dz = 1
    pulselength = 150
    renderstart = 0
    renderstop = 10
     
    # i är tid, j är rum.
    for t in range(1,E.shape[0]): # Iterera över all tid
        if t < pulselength:
            E[t,0] = np.cos(2*np.pi*t*8/200)*np.exp(-(t-pulselength/2)**2 /1e4)
            #E[t,0] = 1
        else :
            E[t,0] = 0
        for z in range(1,E.shape[1]): # För detta tidssteg, iterera över hela rummet
                E[t, z] = E[t-1, z] - dt/dz * (B[t-1,z] - B[t-1, z-1]) 
        for z in range(E.shape[1]-1):
                B[t, z] = B[t-1,z] - dt/dz * (E[t, z+1] - E[t, z])
    for row in range(renderstart, renderstop):
        pyplot.plot(E[row,:], 'b')
        pyplot.ylim(-1, 1)
        pyplot.xlim(0, E.shape[1])
        filedir = 'img/test'
        filename = str(row) + '.png'
        f = filedir+ filename
        pyplot.savefig(f)
        pyplot.clf()


if __name__ == '__main__':
    main()
