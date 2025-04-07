import numpy as np
from parameters import LZ,LR
import scipy.special as sp
from scipy.optimize import fsolve


# Basis function in z using Rational Chebyshev polynomials

def SBz(n, r):
    res = np.sin((2*n+1)*np.arctan(LZ/r)) 
    return res

def dSBz(n,r):
    res = -np.cos((2*n+1)*np.arctan(LZ/r))*(2*n+1)*LZ/(r**2*(1+LZ**2/r**2)) 
    return res


def ddSBz(n,r):
    res = (-np.sin((2*n+1)*np.arctan(LZ/r))*(2*n+1)**2*LZ**2/(r**4*(1+LZ**2/r**2)**2)+
2*np.cos((2*n+1)*np.arctan(LZ/r))*(2*n+1)*LZ/(r**3*(1+LZ**2/r**2))-2*np.cos((2*n+1)*np.arctan(LZ/r))*(2*n+1)*LZ**3/(r**5*(1+LZ**2/r**2)**2))
    return res


# Basis function in r using Rational Chebyshev polynomials

def SBr(n, r):
    res = np.sin((2*n+1)*np.arctan(LR/r)) 
    return res

def dSBr(n,r):
    res = -np.cos((2*n+1)*np.arctan(LR/r))*(2*n+1)*LR/(r**2*(1+LR**2/r**2)) 
    return res


def ddSBr(n,r):
    res = (-np.sin((2*n+1)*np.arctan(LR/r))*(2*n+1)**2*LR**2/(r**4*(1+LR**2/r**2)**2)+
2*np.cos((2*n+1)*np.arctan(LR/r))*(2*n+1)*LR/(r**3*(1+LR**2/r**2))-2*np.cos((2*n+1)*np.arctan(LR/r))*(2*n+1)*LR**3/(r**5*(1+LR**2/r**2)**2))
    return res

