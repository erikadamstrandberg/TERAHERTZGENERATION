#%%
import numpy as np

def Ramp_exp(PLASMASTART,PLASMASTOPP,RAMP_DAMP,Natpunit,Nat,SIZE,dt):
    x = np.array(SIZE)
    i = np.arange(0,SIZE*dt,dt)
    x = 1-np.exp(-i*RAMP_DAMP)
    
    for z in range(SIZE):
        if z > PLASMASTART and z < PLASMASTOPP:
            if z < PLASMASTART + SIZE:
                Nat[z] = x[z-PLASMASTART]*Natpunit
            else:
                Nat[z] = 1*Natpunit
        else:
            Nat[z] = 0
            
def Ramp_linear(RAMPLENGTH,PLASMASTART,PLASMASTOPP,Natpunit,Nat,SIZE,dt):
    x = np.arange(0,RAMPLENGTH*dt,dt)
    xnorm = x/(RAMPLENGTH*dt)
    
    for z in range(SIZE):
        if z > PLASMASTART and z < PLASMASTOPP:
            if z < PLASMASTART + RAMPLENGTH:
                Nat[z] = xnorm[z-PLASMASTART]*Natpunit
            else:
                Nat[z] = 1*Natpunit
        else:
                Nat[z] = 0
