# -*- coding: utf-8 -*-
#

import numpy as np
import matplotlib
matplotlib.use('Agg')
import scipy.constants as const
import matplotlib.pyplot as pyplot

def main():
    maxtime = 500
    maxlength = 500
    dim = [maxtime, maxlength]
    E = np.zeros(dim)
    B = np.zeros(dim)
    dt = 1
    dz = 1
    pulselength = 100
     
    # i är tid, j är rum.
    for t in range(1,len(E)-1): # Iterera över all tid
        for z in range(1,len(E)): # För detta tidssteg, iterera över hela rummet
                E[t, z] = E[t-1, z] - dt/dz * (B[t-1,z] - B[t-1, z-1])
        for z in range(len(B)-1):
                B[t, z] = B[t-1,z] - dt/dz * (E[t, z+1] - E[t, z])
        if t < pulselength:
            E[t,0] = np.cos(np.pi*2*8/200)*np.exp(-(t-pulselength/2)**2/1e4)
            #E[i,0] = 1
        else :
            E[t,0] = 0
    for t in range(100):
        pyplot.plot(E[t,:], 'b')
        pyplot.ylim(-1, 1)
        pyplot.xlim(0, 10)
        filedir = 'img/test'
        filename = str(t) + '.png'
        f = filedir+ filename
        pyplot.savefig(f)

    

if __name__ == '__main__':
    main()
