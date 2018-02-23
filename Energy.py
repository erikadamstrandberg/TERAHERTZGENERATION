# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 15:59:03 2018

@author: Hanna
"""
import scipy.constants as const
import numpy as np

m_e = const.electron_mass
q_e = const.elementary_charge

def Ekin(J,ne):
    Ekin = (m_e/2)*np.abs(J/(ne*q_e))**2
    return Ekin
