import numpy as np
import Plasmaunit as punit

def Micro_i(microREAL,Z,OMEGA_0):
    micro = punit.splasma(microREAL,OMEGA_0)
    micro_cut = np.nonzero(Z > micro)
    micro_i = int(micro_cut[0][0])
    return micro_i