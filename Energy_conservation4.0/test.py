# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 18:55:27 2018

@author: Hanna
"""
#%%
import numpy as np

TIME = 10
SIZE = 10
dt = 1

U_th = np.zeros(TIME)

ne = np.zeros(SIZE)
netemp = np.zeros(SIZE)
kvot2 = np.zeros(SIZE)

for i in range(1,TIME):
    ne = ne+1
    for z in range(SIZE):
        if ne[z]**2 == 0:
            kvot2[z] = 0
        else:
            kvot2[z] = (1/dt)*(ne[z]-netemp[z])   #*(J[z]**2)/(ne[z]**2)
            
    U_th[i] = U_th[i-1] + np.sum(kvot2) 
    netemp = ne