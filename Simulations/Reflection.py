
"""
Created on Thu Feb 22 11:15:31 2018

@author: Hanna
"""
#%%
import scipy.constants as const
import numpy as np
import cmath

q = const.elementary_charge
m = const.electron_mass
epsilon = const.epsilon_0


def R(ne,OMEGA):
    OMEGA_p = np.sqrt(ne*q**2/(m*epsilon))
    EPSILON_r = 1-OMEGA_p**2/OMEGA**2
    n = cmath.sqrt(EPSILON_r)
    R = (np.abs(1-n)/np.abs(1+n))**2
    return R

