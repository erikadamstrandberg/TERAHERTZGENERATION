""" Calculates the propability of a electron ionizing from a atomic core in a external electric field."""
""" The constant names for Landau are choosen to be the same as in Babuskhin.                         """
""" Landau want as arguments a array with the E-field in plasma units at a specific t and the laser   """
""" frequency for the simulation. """

import numpy as np
import scipy.constants as const

#%% Landau

m_e = const.value("electron mass")
e_e = const.value("elementary charge")
epsilon = const.value("electric constant")
hbar = const.value("Planck constant over 2 pi")
eV_to_J = const.value("electron volt-joule relationship")
U_Ar = np.zeros(19)
U_Ar[0] = 0
U_Ar[1] =  15.7596117*eV_to_J # Taget ifrån NIST
U_Ar[2] = 27.62967*eV_to_J
U_Ar[3] = 40.735*eV_to_J
U_Ar[4] = 59.58*eV_to_J
U_Ar[5] = 74.84*eV_to_J
U_Ar[6] = 91.290*eV_to_J
U_Ar[7] = 124.41*eV_to_J
U_Ar[8] = 143.4567*eV_to_J
U_Ar[9] = 422.60*eV_to_J
U_Ar[10] = 479.76*eV_to_J
U_Ar[11] = 540.4*eV_to_J
U_Ar[12] = 619.0*eV_to_J
U_Ar[13] = 685.5*eV_to_J
U_Ar[14] = 755.13*eV_to_J
U_Ar[15] = 855.5*eV_to_J
U_Ar[16] = 918.375*eV_to_J
U_Ar[17] = 4120.6656*eV_to_J
U_Ar[18] = 4426.2228*eV_to_J
U_H = 13.59843449*eV_to_J 

fyra_epsilon_pi = 4*np.pi*epsilon
OMEGA_A = (m_e*e_e**4)/(fyra_epsilon_pi**2*hbar**3)
E_a = (m_e**2*e_e**5)/(fyra_epsilon_pi**3*hbar*4)
r_H = U_Ar/U_H

def Landau(E,LASER_OMEGA,Z):
    W = (4*OMEGA_A*r_H[Z]**(5/2)*(E_a/np.abs(E))*np.exp(-(2*r_H[Z]**(3/2)*(E_a/(3*np.abs(E))))))*LASER_OMEGA
    return W
