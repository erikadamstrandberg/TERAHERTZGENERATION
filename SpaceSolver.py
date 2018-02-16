""" Solves the spatial part of our equation system! """

#%%
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