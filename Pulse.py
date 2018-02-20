""" Create a array with a gaussian pulse. Needed arguements are pulselength and frequency. """
""" Optional are in order amplitude, the mean and the variance of the gaussian envelope    """

import numpy as np

def Gauss(SIZE,OMEGAPRIM,t0,dt):
    El = np.zeros(SIZE)
    START = -SIZE/2
    STOPP = SIZE/2
    t = dt*np.arange(START,STOPP,1)
    for i in range(len(El)):
        El[i] = np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*t0**2)))
    return El

def Envelope(SIZE,t0,dt):
    env = np.zeros(SIZE)
    START = -SIZE/2
    STOPP = SIZE/2
    t = dt*np.arange(START,STOPP,1)
    for i in range(len(env)):
        env[i] = np.exp(-(t[i]**2/(2*t0**2)))
    return env
