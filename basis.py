import time
import numpy as np
import parameters as par
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


np.set_printoptions(precision=16)

t0 = time.process_time()

PR_TEST = par.PR_TEST
PZ_TEST = par.PZ_TEST
LR_TEST = par.LR_TEST
LZ_TEST = par.LZ_TEST


def SBr(n, r):
    res = np.sin((n + 1) * ((np.pi / 2) - np.arctan(r / LR_TEST)))
    return res

def dSBr(n,r):
    res = -np.cos((2*n+1)*np.arctan(LR_TEST/r))*(2*n+1)*LR_TEST/(r**2*(1+LR_TEST**2/r**2)) 
    return res


def ddSBr(n,r):
    res = (-np.sin((2*n+1)*np.arctan(LR_TEST/r))*(2*n+1)**2*LR_TEST**2/(r**4*(1+LR_TEST**2/r**2)**2)+
2*np.cos((2*n+1)*np.arctan(LR_TEST/r))*(2*n+1)*LR_TEST/(r**3*(1+LR_TEST**2/r**2))-2*np.cos((2*n+1)*np.arctan(LR_TEST/r))*(2*n+1)*LR_TEST**3/(r**5*(1+LR_TEST**2/r**2)**2))
    return res

def SBx(n, x):
    r = LR_TEST * x / ( np.sqrt(1 - x**2) )
    res = np.sin((n + 1) * ((np.pi / 2) - np.arctan(r / LR_TEST)))
    return res

def SBz(n, z):
    res = np.sin((n + 1) * ((np.pi / 2) - np.arctan(z / LZ_TEST)))
    return res

def dSBz(n,z):
    res = -np.cos((2*n+1)*np.arctan(LZ_TEST/z))*(2*n+1)*LZ_TEST/(z**2*(1+LZ_TEST**2/z**2)) 
    return res


def ddSBz(n,z):
    res = (-np.sin((2*n+1)*np.arctan(LZ_TEST/z))*(2*n+1)**2*LZ_TEST**2/(z**4*(1+LZ_TEST**2/z**2)**2)+
2*np.cos((2*n+1)*np.arctan(LZ_TEST/z))*(2*n+1)*LZ_TEST/(z**3*(1+LZ_TEST**2/z**2))-2*np.cos((2*n+1)*np.arctan(LZ_TEST/z))*(2*n+1)*LZ_TEST**3/(z**5*(1+LZ_TEST**2/z**2)**2))
    return res

def SBy(n, y):
    z = LZ_TEST * y / ( np.sqrt(1 - y**2) )
    res = np.sin((n + 1) * ((np.pi / 2) - np.arctan(z / LZ_TEST)))
    return res


def SBal(n, r, z):
    list = [SBr(2 * i, r) * SBz(2 * j, z)
            for i in range(PR_TEST + 1) for j in range(PZ_TEST + 1)]
    return list[n]

def drSBal(n, r, z):
    list = [dSBr(2 * i, r) * SBz(2 * j, z)
            for i in range(PR_TEST + 1) for j in range(PZ_TEST + 1)]
    return list[n]

def ddrSBal(n, r, z):
    list = [ddSBr(2 * i, r) * SBz(2 * j, z)
            for i in range(PR_TEST + 1) for j in range(PZ_TEST + 1)]
    return list[n]

def dzSBal(n, r, z):
    list = [SBr(2 * i, r) * dSBz(2 * j, z)
            for i in range(PR_TEST + 1) for j in range(PZ_TEST + 1)]
    return list[n]

def ddzSBal(n, r, z):
    list = [SBr(2 * i, r) * ddSBz(2 * j, z)
            for i in range(PR_TEST + 1) for j in range(PZ_TEST + 1)]
    return list[n]

def SBal_xy(n, x, y):
    list = [SBx(2 * i, x) * SBy(2 * j, y)
            for i in range(PR_TEST + 1) for j in range(PZ_TEST + 1)]
    return list[n]

t1 = time.process_time()
print('Running time (basis):', t1 - t0, 's')