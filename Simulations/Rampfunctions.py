""" Creates a ramp for the atom density to make the current less diskont.    """
""" The exponential ramp scales with the SIZE of the simulation window.      """
""" Linear Erik lets you choose the RAMPLENGTH and adjusts the scaling.      """
""" Linear Hanna lets you choose the RAMP_FRAC for how big of a fraction of  """
""" the atom density to be a ramp.                                           """
  
#%%
import numpy as np

def Ramp_exp(PLASMASTART,PLASMASTOPP,RAMP_DAMP,Natpunit,Nat,SIZE,dz):
    x = np.array(PLASMASTOPP)
    i = np.arange(0,PLASMASTOPP*dz,dz)
    x = 1-np.exp(-i*RAMP_DAMP)
    
    for z in range(SIZE):
        if z > PLASMASTART and z < PLASMASTOPP:
            Nat[z] = x[z-PLASMASTART]*Natpunit
        else:
            Nat[z] = 0
    if SIZE > PLASMASTOPP:
        xback = np.flipud(x)
        for z in range(SIZE):
            if z > PLASMASTART and z < PLASMASTOPP:
                Nat[z] = xback[z]*Nat[z]
            else:
                Nat[z] = 0
        
            
def Ramp_linear_Erik(PLASMASTART,PLASMASTOPP,RAMPLENGTH,Natpunit,Nat,SIZE,dz):
    x = np.arange(0,RAMPLENGTH*dz,dz)
    xnorm = x/(RAMPLENGTH*dz)
    xnormback = np.flipud(xnorm)
    MITTEN = (PLASMASTART + PLASMASTOPP)/2
    
    for z in range(SIZE):
        
        if z > PLASMASTART and z < PLASMASTOPP:  
            if z < PLASMASTART + RAMPLENGTH:
                Nat[z] = xnorm[z-PLASMASTART]*Natpunit
            else:
                Nat[z] = 1*Natpunit   
            if PLASMASTOPP < SIZE:
                if (PLASMASTOPP - RAMPLENGTH) < (PLASMASTART + RAMPLENGTH):
                    if z > MITTEN:
                        Nat[z] = xnormback[z-PLASMASTOPP]*Natpunit
                else:
                    if z > PLASMASTOPP - RAMPLENGTH and z < PLASMASTOPP:
                        Nat[z] = xnormback[z-PLASMASTOPP]*Natpunit
        else:
            Nat[z] = 0
 
           
def Ramp_linear_Hanna(PLASMASTART,PLASMASTOPP,RAMP_FRAC,Natpunit,Nat,SIZE,dz):
    L = PLASMASTOPP-PLASMASTART
    RAMPLENGTH = RAMP_FRAC*L
    x = np.arange(0,RAMPLENGTH*dz,dz)
    xnorm = x/(RAMPLENGTH*dz)
    xnormback = np.flipud(xnorm)
    MITTEN = (PLASMASTART + PLASMASTOPP)/2
    
    for z in range(SIZE):
        
        if z > PLASMASTART and z < PLASMASTOPP:  
            if z < PLASMASTART + RAMPLENGTH:
                Nat[z] = xnorm[z-PLASMASTART]*Natpunit
            else:
                Nat[z] = 1*Natpunit   
            if PLASMASTOPP < SIZE:
                if (PLASMASTOPP - RAMPLENGTH) < (PLASMASTART + RAMPLENGTH):
                    if z > MITTEN:
                        Nat[z] = xnormback[z-PLASMASTOPP]*Natpunit
                else:
                    if z > PLASMASTOPP - RAMPLENGTH and z < PLASMASTOPP:
                        Nat[z] = xnormback[z-PLASMASTOPP]*Natpunit
        else:
            Nat[z] = 0
