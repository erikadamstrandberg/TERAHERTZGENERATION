import scipy.constants as const
import numpy as np

m_e = const.electron_mass
q_e = const.elementary_charge

def Ekin(J,ne):
    Ekin = (m_e/2)*np.abs(J/(ne*q_e))**2
    return Ekin

def Pflux(Espace,Bspace):
    Efs = np.fft.fft(Espace)
    Bfs = np.fft.fft(Bspace)
    Bfstak = np.conj(Bfs)
    
    Sint = np.multiply(Efs,Bfstak)
    Scomplex = (1/2)*sum(Sint)
    Sreal = np.real(Scomplex)
    return Sreal

def EMenergy(E,B):
    W_e = (1/2)*np.sum(E**2)
    W_m = (1/2)*np.sum(B**2)
    W_em = W_e+W_m
    return W_em
