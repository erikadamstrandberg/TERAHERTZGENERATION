#%%

import scipy.constants as const

m_e = const.value("electron mass")
e_e = const.value("elementary charge")
epsilon = const.value("electric constant")
c = const.value("speed of light in vacuum")

def Ereal(E,OMEGA):
    Ereal = E*(m_e*OMEGA*c)/e_e
    return Ereal
    
def Breal(B,OMEGA):
    Breal = B*(m_e*OMEGA)/e_e
    return Breal
    
def Jreal(J,OMEGA):
    Jreal = J*(m_e*OMEGA**2*c*epsilon)/e_e
    return Jreal

def nreal(n,OMEGA):
    nreal = n*(OMEGA**2*m_e*epsilon)/e_e**2
    return nreal