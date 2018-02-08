""" Calculates the propability of a electron ionizing from a atomic core in a external electric field """
""" The constant names for Landau are choosen to be the same as in Babuskhin                          """

import numpy as np
import scipy.constants as const

#%% Landau

m_e = const.value("electron mass")
e_e = const.value("elementary charge")
epsilon = const.value("electric constant")
hbar = const.value("Planck constant over 2 pi")
U_Ar_eV =  15.7596117 # Taget ifrån NIST
U_H_eV = 13.59843449 # Taget ifrån NIST
eV_to_J = const.value("electron volt-joule relationship")
U_Ar_J = U_Ar_eV*eV_to_J
U_H_J = U_H_eV*eV_to_J
ION = U_Ar_J/U_H_J

fyra_epsilon_pi = 4*np.pi*epsilon
OMEGA_A = (m_e*e_e**4)/(fyra_epsilon_pi**2*hbar**3)
E_a = (m_e**2*e_e**5)/(fyra_epsilon_pi**3*hbar*4)
r_H = U_Ar_J/U_H_J

def Landau(E,LASER_OMEGA):
    W = (4*OMEGA_A*r_H**(5/2)*(E_a/np.abs(E))*np.exp(-(2*r_H**(3/2)*(E_a/(3*np.abs(E))))))/LASER_OMEGA
    return W