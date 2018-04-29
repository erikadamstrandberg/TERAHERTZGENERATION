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
    ne[1][:] = ne[0][:]
    return J
 
def N(E,Nat,Ni0,Ni1,Ni2,Ni3,ne,W1,W2,W3,OMEGA_0,dt):
    
    # Calculates all the needed ionization propabilites.
    W1 = Ion.Landau_array_roll(E,OMEGA_0,1)
    W2 = Ion.Landau_array_roll(E,OMEGA_0,2)
    
    # How many neutral atoms are left?
    Ni0[0][:] = Nat-Ni1[0][:]-Ni2[0][:]
    
    # Update eq. for the ionization density.
    Ni1[0][:] = (Ni1[0][:]*(1-(dt/2)*W2)+(dt/2)*W1*(Ni0[0][:]+Ni0[1][:]))/(1+(dt/2)*W2)
    Ni2[0][:] = (Ni2[0][:]+(dt/2)*W2*(Ni1[0][:]+Ni1[1][:]))
    
    # Updates the new electron density.
    ne[0][:] = 1*Ni1[0][:] + 2*Ni2[0][:]
    
    # Makes temp. densities needed for the calculation in the next time step.
    Ni0[1][:] = Ni0[0][:]
    Ni1[1][:] = Ni1[0][:]
    Ni2[1][:] = Ni2[0][:]
    return ne
    
def N4(E,Nat,Ni0,Ni1,Ni2,Ni3,Ni4,Ni5,Ni6,ne,W1,W2,W3,W4,W5,W6,OMEGA_0,dt):
    
    # Calculates all the needed ionization propabilites.
    W1 = Ion.Landau_array_roll(E,OMEGA_0,1)
    W2 = Ion.Landau_array_roll(E,OMEGA_0,2)
    W3 = Ion.Landau_array_roll(E,OMEGA_0,3)
    W4 = Ion.Landau_array_roll(E,OMEGA_0,4)
    W5 = Ion.Landau_array_roll(E,OMEGA_0,5)
    W6 = Ion.Landau_array_roll(E,OMEGA_0,6)
    
    # How many neutral atoms are left?
    Ni0[0][:] = Nat-Ni1[0][:]-Ni2[0][:]-Ni3[0][:]-Ni4[0][:]-Ni5[0][:]-Ni6[0][:]
    
    # Update eq. for the ionization density.
    Ni1[0][:] = (Ni1[0][:]*(1-(dt/2)*W2)+(dt/2)*W1*(Ni0[0][:]+Ni0[1][:]))/(1+(dt/2)*W2)
    Ni2[0][:] = (Ni2[0][:]*(1-(dt/2)*W3)+(dt/2)*W2*(Ni1[0][:]+Ni1[1][:]))/(1+(dt/2)*W3)
    Ni3[0][:] = (Ni3[0][:]*(1-(dt/2)*W4)+(dt/2)*W3*(Ni2[0][:]+Ni2[1][:]))/(1+(dt/2)*W4)
    Ni4[0][:] = (Ni4[0][:]*(1-(dt/2)*W5)+(dt/2)*W4*(Ni3[0][:]+Ni3[1][:]))/(1+(dt/2)*W5)
    Ni5[0][:] = (Ni5[0][:]*(1-(dt/2)*W6)+(dt/2)*W4*(Ni4[0][:]+Ni4[1][:]))/(1+(dt/2)*W6)
    Ni6[0][:] = (Ni6[0][:]+(dt/2)*W6*(Ni5[0][:]+Ni5[1][:]))
    
    # Updates the new electron density.
    ne[0][:] = 1*Ni1[0][:] + 2*Ni2[0][:] + 3*Ni3[0][:] + 4*Ni4[0][:] + 5*Ni5[0][:] + 6*Ni6[0][:]
    
    # Makes temp. densities needed for the calculation in the next time step.
    Ni0[1][:] = Ni0[0][:]
    Ni1[1][:] = Ni1[0][:]
    Ni2[1][:] = Ni2[0][:]
    Ni3[1][:] = Ni3[0][:]
    Ni4[1][:] = Ni4[0][:]
    Ni5[1][:] = Ni5[0][:]
    
    return ne

def N_1e20(E,Nat,Ni0,Ni1,Ni2,Ni3,Ni4,Ni5,Ni6,Ni7,Ni8,ne,W1,W2,W3,W4,W5,W6,W7,W8,OMEGA_0,dt):
    
    # Calculates all the needed ionization propabilites.
    W1 = Ion.Landau_array_roll(E,OMEGA_0,1)
    W2 = Ion.Landau_array_roll(E,OMEGA_0,2)
    W3 = Ion.Landau_array_roll(E,OMEGA_0,3)
    W4 = Ion.Landau_array_roll(E,OMEGA_0,4)
    W5 = Ion.Landau_array_roll(E,OMEGA_0,5)
    W6 = Ion.Landau_array_roll(E,OMEGA_0,6)
    W7 = Ion.Landau_array_roll(E,OMEGA_0,7)
    W8 = Ion.Landau_array_roll(E,OMEGA_0,8)
    
    # How many neutral atoms are left?
    Ni0[0][:] = Nat-Ni1[0][:]-Ni2[0][:]-Ni3[0][:]-Ni4[0][:]-Ni5[0][:]-Ni6[0][:]-Ni7[0][:]-Ni8[0][:]
    
    # Update eq. for the ionization density.
    Ni1[0][:] = (Ni1[0][:]*(1-(dt/2)*W2)+(dt/2)*W1*(Ni0[0][:]+Ni0[1][:]))/(1+(dt/2)*W2)
    Ni2[0][:] = (Ni2[0][:]*(1-(dt/2)*W3)+(dt/2)*W2*(Ni1[0][:]+Ni1[1][:]))/(1+(dt/2)*W3)
    Ni3[0][:] = (Ni3[0][:]*(1-(dt/2)*W4)+(dt/2)*W3*(Ni2[0][:]+Ni2[1][:]))/(1+(dt/2)*W4)
    Ni4[0][:] = (Ni4[0][:]*(1-(dt/2)*W5)+(dt/2)*W4*(Ni3[0][:]+Ni3[1][:]))/(1+(dt/2)*W5)
    Ni5[0][:] = (Ni5[0][:]*(1-(dt/2)*W6)+(dt/2)*W4*(Ni4[0][:]+Ni4[1][:]))/(1+(dt/2)*W6)
    Ni6[0][:] = (Ni6[0][:]*(1-(dt/2)*W7)+(dt/2)*W6*(Ni5[0][:]+Ni5[1][:]))/(1+(dt/2)*W7)
    Ni7[0][:] = (Ni7[0][:]*(1-(dt/2)*W8)+(dt/2)*W7*(Ni6[0][:]+Ni6[1][:]))/(1+(dt/2)*W8)
    Ni8[0][:] = (Ni8[0][:]+(dt/2)*W8*(Ni7[0][:]+Ni7[1][:]))
    
    # Updates the new electron density.
    ne[0][:] = 1*Ni1[0][:] + 2*Ni2[0][:] + 3*Ni3[0][:] + 4*Ni4[0][:] + 5*Ni5[0][:] + 6*Ni6[0][:] + 7*Ni7[0][:] + 8*Ni8[0][:]
    
    # Makes temp. densities needed for the calculation in the next time step.
    Ni0[1][:] = Ni0[0][:]
    Ni1[1][:] = Ni1[0][:]
    Ni2[1][:] = Ni2[0][:]
    Ni3[1][:] = Ni3[0][:]
    Ni4[1][:] = Ni4[0][:]
    Ni5[1][:] = Ni5[0][:]
    Ni6[1][:] = Ni6[0][:]
    Ni7[1][:] = Ni7[0][:]
    Ni8[1][:] = Ni8[0][:]
    return ne
    
def N_full(E,Nat,Ni0,Ni1,Ni2,Ni3,Ni4,Ni5,Ni6,Ni7,Ni8,Ni9,Ni10,Ni11,Ni12,Ni13,Ni14,Ni15,Ni16,Ni17,Ni18,ne,W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,W14,W15,W16,W17,W18,OMEGA_0,dt):
    
    # Calculates all the needed ionization propabilites.
    W1 = Ion.Landau_array_roll(E,OMEGA_0,1)
    W2 = Ion.Landau_array_roll(E,OMEGA_0,2)
    W3 = Ion.Landau_array_roll(E,OMEGA_0,3)
    W4 = Ion.Landau_array_roll(E,OMEGA_0,4)
    W5 = Ion.Landau_array_roll(E,OMEGA_0,5)
    W6 = Ion.Landau_array_roll(E,OMEGA_0,6)
    W7 = Ion.Landau_array_roll(E,OMEGA_0,7)
    W8 = Ion.Landau_array_roll(E,OMEGA_0,8)
    W9 = Ion.Landau_array_roll(E,OMEGA_0,9)
    W10 = Ion.Landau_array_roll(E,OMEGA_0,10)
    W11 = Ion.Landau_array_roll(E,OMEGA_0,11)
    W12 = Ion.Landau_array_roll(E,OMEGA_0,12)
    W13 = Ion.Landau_array_roll(E,OMEGA_0,13)
    W14 = Ion.Landau_array_roll(E,OMEGA_0,14)
    W15 = Ion.Landau_array_roll(E,OMEGA_0,15)
    W16 = Ion.Landau_array_roll(E,OMEGA_0,16)
    W17 = Ion.Landau_array_roll(E,OMEGA_0,17)
    W18 = Ion.Landau_array_roll(E,OMEGA_0,18)
    
    # How many neutral atoms are left?
    Ni0[0][:] = Nat-Ni1[0][:]-Ni2[0][:]-Ni3[0][:]-Ni4[0][:]-Ni5[0][:]-Ni6[0][:]-Ni7[0][:]-Ni8[0][:]-Ni9[0][:]-Ni10[0][:]-Ni11[0][:]-Ni12[0][:]-Ni13[0][:]-Ni14[0][:]-Ni15[0][:]-Ni16[0][:]-Ni17[0][:]-Ni18[0][:]
    
    # Update eq. for the ionization density.
    Ni1[0][:] = (Ni1[0][:]*(1-(dt/2)*W2)+(dt/2)*W1*(Ni0[0][:]+Ni0[1][:]))/(1+(dt/2)*W2)    
    Ni2[0][:] = (Ni2[0][:]*(1-(dt/2)*W3)+(dt/2)*W2*(Ni1[0][:]+Ni1[1][:]))/(1+(dt/2)*W3)
    Ni3[0][:] = (Ni3[0][:]*(1-(dt/2)*W4)+(dt/2)*W3*(Ni2[0][:]+Ni2[1][:]))/(1+(dt/2)*W4)
    Ni4[0][:] = (Ni4[0][:]*(1-(dt/2)*W5)+(dt/2)*W4*(Ni3[0][:]+Ni3[1][:]))/(1+(dt/2)*W5)
    Ni5[0][:] = (Ni5[0][:]*(1-(dt/2)*W6)+(dt/2)*W5*(Ni4[0][:]+Ni4[1][:]))/(1+(dt/2)*W6)
    Ni6[0][:] = (Ni6[0][:]*(1-(dt/2)*W7)+(dt/2)*W6*(Ni5[0][:]+Ni5[1][:]))/(1+(dt/2)*W7)
    Ni7[0][:] = (Ni7[0][:]*(1-(dt/2)*W8)+(dt/2)*W7*(Ni6[0][:]+Ni6[1][:]))/(1+(dt/2)*W8)
    Ni8[0][:] = (Ni8[0][:]*(1-(dt/2)*W9)+(dt/2)*W8*(Ni7[0][:]+Ni7[1][:]))/(1+(dt/2)*W9)
    Ni9[0][:] = (Ni9[0][:]*(1-(dt/2)*W10)+(dt/2)*W9*(Ni8[0][:]+Ni8[1][:]))/(1+(dt/2)*W10)
    Ni10[0][:] = (Ni10[0][:]*(1-(dt/2)*W11)+(dt/2)*W10*(Ni9[0][:]+Ni9[1][:]))/(1+(dt/2)*W11)
    Ni11[0][:] = (Ni11[0][:]*(1-(dt/2)*W12)+(dt/2)*W11*(Ni10[0][:]+Ni10[1][:]))/(1+(dt/2)*W12)
    Ni12[0][:] = (Ni12[0][:]*(1-(dt/2)*W13)+(dt/2)*W12*(Ni11[0][:]+Ni11[1][:]))/(1+(dt/2)*W13)
    Ni13[0][:] = (Ni13[0][:]*(1-(dt/2)*W14)+(dt/2)*W13*(Ni12[0][:]+Ni12[1][:]))/(1+(dt/2)*W14)
    Ni14[0][:] = (Ni14[0][:]*(1-(dt/2)*W15)+(dt/2)*W14*(Ni13[0][:]+Ni13[1][:]))/(1+(dt/2)*W15)
    Ni15[0][:] = (Ni15[0][:]*(1-(dt/2)*W16)+(dt/2)*W15*(Ni14[0][:]+Ni14[1][:]))/(1+(dt/2)*W16)
    Ni16[0][:] = (Ni16[0][:]*(1-(dt/2)*W17)+(dt/2)*W16*(Ni15[0][:]+Ni15[1][:]))/(1+(dt/2)*W17)
    Ni17[0][:] = (Ni17[0][:]*(1-(dt/2)*W18)+(dt/2)*W17*(Ni16[0][:]+Ni16[1][:]))/(1+(dt/2)*W18)
    Ni18[0][:] = (Ni18[0][:]+(dt/2)*W18*(Ni17[0][:]+Ni17[1][:]))
    
    # Updates the new electron density.
    ne[0][:] = 1*Ni1[0][:] + 2*Ni2[0][:] + 3*Ni3[0][:] + 4*Ni4[0][:] + 5*Ni5[0][:] + 6*Ni6[0][:] + 7*Ni7[0][:] + 8*Ni8[0][:] + 9*Ni9[0][:] + 10*Ni10[0][:] + 11*Ni11[0][:] + 12*Ni12[0][:] + 13*Ni13[0][:] + 14*Ni14[0][:] + 15*Ni15[0][:] + 16*Ni16[0][:] + 17*Ni17[0][:] + 18*Ni18[0][:]
    
    # Makes temp. densities needed for the calculation in the next time step.
    Ni0[1][:] = Ni0[0][:]
    Ni1[1][:] = Ni1[0][:]
    Ni2[1][:] = Ni2[0][:]
    Ni3[1][:] = Ni3[0][:]
    Ni4[1][:] = Ni4[0][:]
    Ni5[1][:] = Ni5[0][:]
    Ni6[1][:] = Ni6[0][:]
    Ni7[1][:] = Ni7[0][:]
    Ni8[1][:] = Ni8[0][:]
    Ni9[1][:] = Ni9[0][:]
    Ni10[1][:] = Ni10[0][:]
    Ni11[1][:] = Ni11[0][:]
    Ni12[1][:] = Ni12[0][:]
    Ni13[1][:] = Ni13[0][:]
    Ni14[1][:] = Ni14[0][:]
    Ni15[1][:] = Ni15[0][:]
    Ni16[1][:] = Ni16[0][:]
    Ni17[1][:] = Ni17[0][:]
    Ni18[1][:] = Ni18[0][:]
    return ne
