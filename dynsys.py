import numpy as np
import collocations as col
import parameters as par
import time
t0 = time.process_time()

np.set_printoptions(precision=16)

PZ = par.PZ_TEST
PR = par.PR_TEST
# ddzSB_cyl = col.ddzSB
# drSB = col.drSB
# ddrSB_cyl = col.ddrSB
col_r_shift = col.col_r_shift
SBal = col.SBal
SBal_inv = col.SBal_inv
ddzSBal = col.ddzSBal
ddrSBal = col.ddrSBal
drSBal = col.drSBal

r = np.tile(col_r_shift.repeat(PZ + 1).reshape(-1, 1), (1, (PR + 1) * (PZ + 1)))

def dda(c):

    # drphi = np.dot(drSBphi,c)
    # drrphi = np.dot(ddrSBphi,c)
    # dzzphi = np.dot(ddzSBphi,c)

    # RHS = drrphi + drphi*1/r  + dzzphi

    # return np.dot(RHS,SBphi_inv)
    return np.dot(np.dot(c, ddrSBal + drSBal*1/r  + ddzSBal),SBal_inv)

def phi(c):
    return np.dot(SBal,c)

t1 = time.process_time()
print('Running time (dynamical system):', t1 - t0, 's')

