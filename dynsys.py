import numpy as np
import collocations as col
import parameters as par

PZ = par.PZ
ddzSB_cyl = col.ddzSB
drSB = col.drSB
ddrSB_cyl = col.ddrSB
rcol = col.rcol
SB = col.SB
SB_inv = col.SB_inv
SB_inv = col.SB_inv
ddzSB = col.ddzSB
drSB = col.drSB
ddrSB = col.ddrSB

rrep = np.repeat(rcol,PZ+1)

def dda(c):

    drphi = np.dot(c,drSB)
    drrphi = np.dot(c,ddrSB)
    dzzphi = np.dot(c,ddzSB)

    RHS = drrphi + drphi*1/rrep  + dzzphi

    return np.dot(RHS,SB_inv)

def phi(c):
    return np.dot(c,SB)

