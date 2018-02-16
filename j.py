# -*- coding: utf-8 -*-

import numpy as np

def J(E, n, v, dt, dz, jstart, jstartlength):
    "Returns a current field J(t, z) of same dimension as E."
    # E field, n field, v int, dt int, dz int, jstart start value, jstartlength time the start value is to be set
    J = np.array(np.shape(E))
    for t in range(1, J.shape[0]-1):
        if t < jstartlength:
            J[t, 0] = jstart
        for z in range(J.shape[1]):
            J[t, z] = (1 - v*dt/2)/(1 + v*dt/2)*J[t-1, z] + dt/2*(n[t-1, z] + n[t, z])*E[t-1, z]/(1 + v*dt/2)
                
    return J        
