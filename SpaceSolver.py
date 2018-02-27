""" Solves the spatial part of our equation system! """

#%%
import Ionization as Ion

def E(E,B,J,dt,dz):
    for z in range(1,len(E)):
        E[z] = E[z]-(dt/dz)*(B[z]-B[z-1])-dt*J[z-1]
    return E
    
def B(E,B,dt,dz):
    for z in range(len(B)-1):
        B[z] = B[z]-(dt/dz)*(E[z+1]-E[z])
    return B

def J(E,J,n,ntemp,nu,dt,dz):
    for z in range(len(J)):
        J[z] = ((1 - nu*dt/2)*J[z]+dt/2*(n[z]+ntemp[z])*E[z])/(1+nu*dt/2)
    return J

def N(E,Nat,Ni0,Ni1,Ni2,Ni0temp,Ni1temp,ne,W1,W2,W3,OMEGA_0,dt):
    for z in range(len(E)):
        W1[z] = Ion.Landau_element(E[z],OMEGA_0,1,dt)
        W2[z] = Ion.Landau_element(E[z],OMEGA_0,2,dt)
        W3[z] = Ion.Landau_element(E[z],OMEGA_0,3,dt)
        Ni0[z] = Nat[z]-Ni1[z]-Ni2[z]
        Ni1[z] = (Ni1[z]*(1-(dt/2)*W2[z])+(dt/2)*W1[z]*(Ni0[z]+Ni0temp[z]))/(1+(dt/2)*W2[z])
        Ni2[z] = (Ni2[z]*(1-(dt/2)*W3[z])+(dt/2)*W2[z]*(Ni1[z]+Ni1temp[z]))/(1+(dt/2)*W3[z])
        ne[z] = 1*Ni1[z] + 2*Ni2[z]
return ne
