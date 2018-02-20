""" Create an array with a gaussian pulse. Needed arguements are pulselength and frequency. """
""" Optional are in order amplitude, the mean and the variance of the gaussian envelope    """

import numpy as np
    
def Gauss_forward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,VAR,dt):
    El = np.zeros(PULSELENGTH)
    Bl = np.zeros(PULSELENGTH)
    START = -PULSELENGTH/2
    STOPP = PULSELENGTH/2
    t = dt*np.arange(START,STOPP,1)
    for i in range(len(El)):
        El[i] = E0*np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*VAR)))/np.sqrt(2*np.pi*VAR)
    for i in range(len(Bl)):
        Bl[i] = E0*(np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*VAR))))/np.sqrt(2*np.pi*VAR)
    E[PULSESTART:PULSESTART+PULSELENGTH] = El
    B[PULSESTART:PULSESTART+PULSELENGTH] = Bl
    
def Gauss_backward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,VAR,dt):
    El = np.zeros(PULSELENGTH)
    Bl = np.zeros(PULSELENGTH)
    START = -PULSELENGTH/2
    STOPP = PULSELENGTH/2
    t = dt*np.arange(START,STOPP,1)
    for i in range(len(El)):
        El[i] = E0*np.cos(OMEGAPRIM*t[i-1])*np.exp(-(t[i-1]**2/(2*VAR)))/np.sqrt(2*np.pi*VAR)
    for i in range(len(Bl)):
        Bl[i] = -E0*(np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*VAR))))/np.sqrt(2*np.pi*VAR)
    E[PULSESTART:PULSESTART+PULSELENGTH] = El
    B[PULSESTART:PULSESTART+PULSELENGTH] = Bl
