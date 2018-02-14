""" Converts your value from plasmaunits and back! """
""" Needs the OMEGA from your simulation           """

#%%
import scipy.constants as const

m_e = const.value("electron mass")
e_e = const.value("elementary charge")
epsilon = const.value("electric constant")
c = const.value("speed of light in vacuum")

def Ereal(E,OMEGA):
    Ereal = E*(m_e*OMEGA*c)/e_e
    return Ereal

def Eplasma(E,OMEGA):
    Eplasma = E*(e_e/(m_e*OMEGA*c))
    return Eplasma
    
def Breal(B,OMEGA):
    Breal = B*(m_e*OMEGA)/e_e
    return Breal

def Bplasma(B,OMEGA):
    Bplasma = B*(e_e/(m_e*OMEGA))
    return Bplasma
    
def Jreal(J,OMEGA):
    Jreal = J*(m_e*OMEGA**2*c*epsilon)/e_e
    return Jreal

def Jplasma(J,OMEGA):
    Jplasma = J*(e_e/(m_e*OMEGA**2*c*epsilon))
    return Jplasma
    
def nreal(n,OMEGA):
    nreal = n*(OMEGA**2*m_e*epsilon)/e_e**2
    return nreal

def nplasma(n,OMEGA):
    nplasma = n*(e_e**2/(OMEGA**2*m_e*epsilon))
    return nplasma

def treal(t,OMEGA):
    treal = t/OMEGA
    return treal
    
def tplasma(t,OMEGA):
    tplasma = t*OMEGA
    return tplasma

def sreal(s,OMEGA):
    sreal = s*(c/OMEGA)
    return sreal
    
def splasma(s,OMEGA):
    splasma = s*(OMEGA/c)
    return splasma

def omegareal(OMEGAPRIM,OMEGA):
    omegareal = OMEGAPRIM*OMEGA
    return omegareal
    
def omegaplasma(OMEGAPRIM,OMEGA):
    omegaplasma = OMEGAPRIM/OMEGA
    return omegaplasma
