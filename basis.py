import numpy as np
from parameters import LX,LR,PR,PX
import scipy.special as sp

# Basis function in r using Odd Chebyshev Polynomials

def SB(n, r):
    res = np.sin((n+1)*np.arctan(LR/r)) 
    return res

def dSB(n,r):
    res = -np.cos((n+1)*np.arctan(LR/r))*(n+1)*LR/(r**2*(1+LR**2/r**2)) 
    return res


def ddSB(n,r):
    res = (-np.sin((n+1)*np.arctan(LR/r))*(n+1)**2*LR**2/(r**4*(1+LR**2/r**2)**2)+
2*np.cos((n+1)*np.arctan(LR/r))*(n+1)*LR/(r**3*(1+LR**2/r**2))-2*np.cos((n+1)*np.arctan(LR/r))*(n+1)*LR**3/(r**5*(1+LR**2/r**2)**2))
    return res

# Basis function in x using Legendre Polynomials

def P(i, x):
    return sp.legendre(i)(x)

def dP(i, x):
    return sp.legendre(i).deriv()(x)

def ddP(i, x):
    return sp.legendre(i).deriv().deriv()(x)

###Big Basis

def Basis(n, r, x):
    list = [SB(2 * i, r) * P(2 * j, x)
            for i in range(PR + 1) for j in range(PX + 1)]
    return list[n]

# def Basis(i, j, r, x):
#     list = SB(2 * i, r) * P(2 * j, x)
#     return list