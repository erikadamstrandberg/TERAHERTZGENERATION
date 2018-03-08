""" initializes the E and B field for a only forward or backward propagation pulse! """
""" NOTE: Need the E-field to be calculated first!                                  """

import numpy as np
    
def Gauss_forward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt):
    El = np.zeros(PULSELENGTH)
    Bl = np.zeros(PULSELENGTH)
    START = -PULSELENGTH/2  
    STOPP = PULSELENGTH/2
    t = dt*np.arange(START,STOPP,1)
    for i in range(len(El)):
        El[i] = E0*np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*t0**2)))
    for i in range(len(Bl)):
        Bl[i] = E0*(np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*t0**2))))
    E[PULSESTART:PULSESTART+PULSELENGTH] = El
    B[PULSESTART:PULSESTART+PULSELENGTH] = Bl
    
def Gauss_backward(E,B,E0,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt):
    El = np.zeros(PULSELENGTH)
    Bl = np.zeros(PULSELENGTH)
    START = -PULSELENGTH/2
    STOPP = PULSELENGTH/2
    t = dt*np.arange(START,STOPP,1)
    for i in range(len(El)):
        El[i] = E0*np.cos(OMEGAPRIM*t[i-1])*np.exp(-(t[i-1]**2/(2*t0**2)))
    for i in range(len(Bl)):
        Bl[i] = -E0*(np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*t0**2))))
    E[PULSESTART:PULSESTART+PULSELENGTH] = El
    B[PULSESTART:PULSESTART+PULSELENGTH] = Bl
   
def TwoC_forward(E,B,E0,xi,phi,PULSELENGTH,PULSESTART,OMEGAPRIM,t0,dt):
    El = np.zeros(PULSELENGTH)
    Bl = np.zeros(PULSELENGTH)
    START = -PULSELENGTH/2  
    STOPP = PULSELENGTH/2
    t = dt*np.arange(START,STOPP,1)
    for i in range(len(El)):
        El[i] = E0*(np.sqrt(1-xi)*np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*t0**2)))+np.sqrt(xi)*np.cos(2*OMEGAPRIM*t[i]+phi)*np.exp(-4*(t[i]**2/(2*t0**2))))
    for i in range(len(Bl)):
        Bl[i] = E0*(np.sqrt(1-xi)*np.cos(OMEGAPRIM*t[i])*np.exp(-(t[i]**2/(2*t0**2)))+np.sqrt(xi)*np.cos(2*OMEGAPRIM*t[i]+phi)*np.exp(-4*(t[i]**2/(2*t0**2))))
    E[PULSESTART:PULSESTART+PULSELENGTH] = El
    B[PULSESTART:PULSESTART+PULSELENGTH] = Bl
