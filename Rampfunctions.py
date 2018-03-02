""" Creates a ramp for the atom density to make the current less diskont.    """
""" The exponential ramp scales with the SIZE of the simulation window       """
""" With the linear ramp you have to scale the RAMPLENGTH with the SIZE      """

#%%
import numpy as np

def Ramp_exp(PLASMASTART,PLASMASTOPP,RAMP_DAMP,Natpunit,Nat,SIZE,dz):
    x = np.array(SIZE)
    i = np.arange(0,SIZE*dz,dz)
    x = 1-np.exp(-i*RAMP_DAMP)
    
    for z in range(SIZE):
        if z > PLASMASTART and z < PLASMASTOPP:
            if z < PLASMASTART + SIZE:
                Nat[z] = x[z-PLASMASTART]*Natpunit
            else:
                Nat[z] = 1*Natpunit
        else:
            Nat[z] = 0
            
def Ramp_linear(RAMPLENGTH,PLASMASTART,PLASMASTOPP,Natpunit,Nat,SIZE,dz):
    x = np.arange(0,RAMPLENGTH*dZ,dz)
    xnorm = x/(RAMPLENGTH*dz)
    
    for z in range(SIZE):
        if z > PLASMASTART and z < PLASMASTOPP:
            if z < PLASMASTART + RAMPLENGTH:
                Nat[z] = xnorm[z-PLASMASTART]*Natpunit
            else:
                Nat[z] = 1*Natpunit
        else:
                Nat[z] = 0
