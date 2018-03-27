""" Solves the spatial part of our equation system! """
""" Remember to solve in order E,B,ne,J             """

#%%
import Ionization as Ion
import numpy as np

def E(E,B,J,dt,dz):
    E = E-(dt/dz)*(B-np.roll(B,1))-dt*J
    return E
    
def B(E,B,dt,dz):
    B = B-(dt/dz)*(np.roll(E,-1)-E)
    return B

def J(E,J,ne,nu,dt,dz):
        
    # Update eq. for the J-field.
    J = ((1-nu*dt/2)*J+(dt/2)*(ne[0]+ne[1])*E)/(1+nu*dt/2)
    
    # Makes temp. electron density need for the next time step.
    ne[1] = ne[0]
    return J

def N(E,Nat,Ni0,Ni1,Ni2,ne,W1,W2,W3,OMEGA_0,dt):
    
    # Calculates all the needed ionization propabilites.
    W1 = Ion.Landau_array(E,W1,OMEGA_0,1,dt)
    W2 = Ion.Landau_array(E,W2,OMEGA_0,2,dt)
    W3 = Ion.Landau_array(E,W2,OMEGA_0,3,dt)
    
    # How many neutral atoms are left?
    Ni0[0] = Nat-Ni1[0]-Ni2[0]
    
    # Update eq. for the ionization density.
    Ni1[0] = (Ni1[0]*(1-(dt/2)*W2)+(dt/2)*W1*(Ni0[0]+Ni0[1]))/(1+(dt/2)*W2)
    Ni2[0] = (Ni2[0]*(1-(dt/2)*W3)+(dt/2)*W2*(Ni1[0]+Ni1[1]))/(1+(dt/2)*W3)
    
    # Updates the new electron density.
    ne[0] = 1*Ni1[0] + 2*Ni2[0]
    
    # Makes temp. densities needed for the calculation in the next time step.
    Ni0[1] = Ni0[0]
    Ni1[1] = Ni1[0]
    Ni2[1] = Ni2[0]
    return ne
