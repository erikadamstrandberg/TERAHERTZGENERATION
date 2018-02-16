""" Create a array with a gaussian pulse. Needed arguements are pulselength and frequency. """
""" Optional are in order amplitude, the mean and the variance of the gaussian envelope    """

import numpy as np
def Laser_gauss_cos(SIZE,OMEGA,*OPTIONS):
    El = np.zeros(SIZE)
    if not OPTIONS:
        for i in range(len(El)):
            El[i] = np.cos(OMEGA*i)*np.exp(-(i-SIZE/2)**2/1e4)
    elif len(OPTIONS) == 1:
        AMPLITUDE = OPTIONS[0]
        for i in range(len(El)):
            El[i] = AMPLITUDE*np.cos(OMEGA*i)*np.exp(-(i-SIZE/2)**2/1e4)
    elif len(OPTIONS) == 2:
        AMPLITUDE = OPTIONS[0]
        MEAN = OPTIONS[1]
        for i in range(len(El)):
            El[i] = AMPLITUDE*np.cos(OMEGA*i)*np.exp(-(i-MEAN)**2/1e4)
    elif len(OPTIONS) == 3:
        AMPLITUDE = OPTIONS[0]
        MEAN = OPTIONS[1]
        VAR = OPTIONS[2]
        for i in range(len(El)):
            El[i] = AMPLITUDE*np.cos(OMEGA*i)*np.exp(-(i-MEAN)**2/VAR)
    else:
        for i in range(len(El)):
            El[i] = np.cos(OMEGA*i)*np.exp(-(i-SIZE/2)**2/1e4)

    return El
    
def Laser_gauss_sin(SIZE,OMEGA,*OPTIONS):
    El = np.zeros(SIZE)
    if not OPTIONS:
        for i in range(len(El)):
            El[i] = np.sin(OMEGA*i)*np.exp(-(i-SIZE/2)**2/1e4)
    elif len(OPTIONS) == 1:
        MEAN = OPTIONS[0]
        for i in range(len(El)):
            El[i] = np.sin(OMEGA*i)*np.exp(-(i-MEAN)**2/1e4)
    elif len(OPTIONS) == 2:
        MEAN = OPTIONS[0]
        VAR = OPTIONS[1]
        for i in range(len(El)):
            El[i] = np.sin(OMEGA*i)*np.exp(-(i-MEAN)**2/VAR)
    else:
        for i in range(len(El)):
            El[i] = np.sin(OMEGA*i)*np.exp(-(i-SIZE/2)**2/1e4)

    return El