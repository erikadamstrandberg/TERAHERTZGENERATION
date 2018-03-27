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
    return J + dt*ne[0,:]*E
    #dt2 = dt/2
    #for z in range(len(J)):
        
        # Update eq. for the J-field.
        #J[z] = ((1-nu*dt2)*J[z]+(dt2)*(ne[0][z]+ne[1][z])*E[z])/(1+nu*dt2)
        
        # Makes temp. electron density need for the next time step.
        #ne[1][z] = ne[0][z]

def N(E,Nat,Ni0,Ni1,Ni2,ne,W1,W2,W3,OMEGA_0,dt):
    
    dt2 = dt/2
    for z in range(len(E)):
        
        # Calculates all the needed ionization propabilites.
        W1[z] = Ion.Landau_element(E[z],OMEGA_0,1,dt)
        W2[z] = Ion.Landau_element(E[z],OMEGA_0,2,dt)
        W3[z] = Ion.Landau_element(E[z],OMEGA_0,3,dt)
        
        # How many neutral atoms are left?
        Ni0[0][z] = Nat[z]-Ni1[0][z]-Ni2[0][z]
        
        # Update eq. for the ionization density.
        Ni1[0][z] = (Ni1[0][z]*(1-(dt2)*W2[z])+(dt2)*W1[z]*(Ni0[0][z]+Ni0[1][z]))/(1+(dt2)*W2[z])
        Ni2[0][z] = (Ni2[0][z]*(1-(dt2)*W3[z])+(dt2)*W2[z]*(Ni1[0][z]+Ni1[1][z]))/(1+(dt2)*W3[z])
        
        # Updates the new electron density.
        ne[0][z] = 1*Ni1[0][z] + 2*Ni2[0][z]
        
        # Makes temp. densities needed for the calculation in the next time step.
        Ni0[1][z] = Ni0[0][z]
        Ni1[1][z] = Ni1[0][z]
        Ni2[1][z] = Ni2[0][z]
    return ne
