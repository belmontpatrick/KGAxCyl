import numpy as np
import collocation as cp
import collocation_matrices as cm

rcol = cp.rcol
xcol = cp.xcol
Psi = cm.Psi
Psi_inv = cm.Psi_inv
drPsi = cm.drPsi
ddrPsi = cm.ddrPsi
dxPsi = cm.dxPsi
ddxPsi = cm.ddxPsi

Rm = np.repeat(rcol,len(xcol))
Rr = Rm.reshape(-1,1)

Xm = np.repeat(xcol,len(rcol))
Xr = Xm.reshape(-1,1)

def dda(c):

    drphi = np.dot(drPsi,c)
    drrphi = np.dot(ddrPsi,c)
    dxxphi = np.dot(ddxPsi,c)
    dxphi = np.dot(dxPsi,c)


    RHS = drrphi + drphi*2/Rr  + (1/ ( Rr ) **2) * ((1 - Xr **2) * dxxphi * 2 * Xr * dxphi)

    return np.dot(Psi_inv,RHS)

def phi(c):
    return np.dot(Psi,c)
