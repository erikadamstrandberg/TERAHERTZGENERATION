""" initializes the E and B field for a only forward or backward propagation pulse! """
""" NOTE: Need the E-field to be calculated first!                                  """
""" NOTE: PULSELENGTH and PULSSTART in terms of 'index units' """

import numpy as np
    
def Gauss_forward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt):
    El = np.zeros(PULSELENGTH+1)
    Bl = np.zeros(PULSELENGTH+1)
    START = -PULSELENGTH/2  
    STOPP = PULSELENGTH/2+1
    t = dt*np.arange(START,STOPP,1)
    for i in range(len(El)):
        El[i] = E0*np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(t0**2)))
    for i in range(len(Bl)):
        Bl[i] = E0*(np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(t0**2))))
    E[PULSESTART:PULSESTART+PULSELENGTH+1] = El
    B[PULSESTART:PULSESTART+PULSELENGTH+1] = Bl
    
def Gauss_backward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt):
    El = np.zeros(PULSELENGTH+1)
    Bl = np.zeros(PULSELENGTH+1)
    START = -PULSELENGTH/2
    STOPP = PULSELENGTH/2+1
    t = dt*np.arange(START,STOPP,1)
    for i in range(len(El)):
        El[i] = E0*np.cos(OMEGAPRIM*t[i-1])*np.exp(-(t[i-1]**2/(t0**2)))
    for i in range(len(Bl)):
        Bl[i] = -E0*(np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(t0**2))))
    E[PULSESTART:PULSESTART+PULSELENGTH+1] = El
    B[PULSESTART:PULSESTART+PULSELENGTH+1] = Bl
