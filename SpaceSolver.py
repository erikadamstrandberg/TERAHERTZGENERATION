""" Solves the spatial part of our equation system! """
""" Remember to solve in order E,B,ne,J             """

#%%
import Ionization as Ion

def E(E,B,J,dt,dz):
    for z in range(1,len(E)):
        
        # Update eq. for the E-field.
        E[z] = E[z]-(dt/dz)*(B[z]-B[z-1])-dt*J[z]
    return E
    
def B(E,B,dt,dz):
    for z in range(len(B)-1):
        
        # Update eq. for the B-field.
        B[z] = B[z]-(dt/dz)*(E[z+1]-E[z])
    return B

def J(E,J,ne,nu,dt,dz):
    for z in range(len(J)-1):
        
        # Update eq. for the J-field.
        J[z] = ((1-nu*dt/2)*J[z]+(dt/2)*(ne[0][z]+ne[1][z])*E[z+1])/(1+nu*dt/2)
        
        # Makes temp. electron density need for the next time step.
        ne[1][z] = ne[0][z]
    return J

def N(E,Nat,Ni0,Ni1,Ni2,ne,W1,W2,W3,OMEGA_0,dt):
    for z in range(len(E)-1):
        
        # Calculates all the needed ionization propabilites.
        W1[z] = Ion.Landau_element(E[z+1],OMEGA_0,1,dt)
        W2[z] = Ion.Landau_element(E[z+1],OMEGA_0,2,dt)
        W3[z] = Ion.Landau_element(E[z+1],OMEGA_0,3,dt)
        
        # How many neutral atoms are left?
        Ni0[0][z] = Nat[z]-Ni1[0][z]-Ni2[0][z]
        
        # Update eq. for the ionization density.
        Ni1[0][z] = (Ni1[0][z]*(1-(dt/2)*W2[z])+(dt/2)*W1[z]*(Ni0[0][z]+Ni0[1][z]))/(1+(dt/2)*W2[z])
        Ni2[0][z] = (Ni2[0][z]*(1-(dt/2)*W3[z])+(dt/2)*W2[z]*(Ni1[0][z]+Ni1[1][z]))/(1+(dt/2)*W3[z])
        
        # Updates the new electron density.
        ne[0][z] = 1*Ni1[0][z] + 2*Ni2[0][z]
        
        # Makes temp. densities needed for the calculation in the next time step.
        Ni0[1][z] = Ni0[0][z]
        Ni1[1][z] = Ni1[0][z]
        Ni2[1][z] = Ni2[0][z]
    return ne
