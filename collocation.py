import numpy as np
import basis as bs
import parameters as par
import scipy.special as sp
from scipy.optimize import fsolve

PR = par.PR
LR = par.LR
PX = par.PX
P = bs.P
dP = bs.dP
ddP = bs.ddP
SB = bs.SB
dSB = bs.dSB
ddSB = bs.ddSB


"Legendre Collocation points in x"

def P_colpoints(n):
    P_prime = sp.legendre(2 * PX + 3).deriv()
    x_roots = fsolve(P_prime, np.cos(np.pi * (np.arange(1, 2 * PX + 3) / (2 * PX + 3))))
    # x_col_prel = np.sort(x_roots)
    x_col = x_roots[:PX + 1]
    # x_col = np.flip(x_col_prel[:PX + 1])

    return x_col

xcol = P_colpoints(PX)

# print("xcol", xcol)


"Base Matrix (Legendre Polinomials) in x: "

P_x = np.zeros([PX+1,PX+1])
dP_x = np.zeros([PX+1,PX+1])
ddP_x = np.zeros([PX+1,PX+1])

for i in range(PX+1):
  P_x[i,] = P(2*i,xcol)

for i in range(PX+1):
  dP_x[i,] = dP(2*i,xcol)

for i in range(PX+1):
  ddP_x[i,] = ddP(2*i,xcol)

P_x_inv = np.linalg.inv(P_x)

#print(xcol)

#np.savetxt("xcolPX2.txt", xcol, fmt="%.16e", delimiter="\t")

# print("xcol_Henrique", xcol_Henrique)
def SB_colpoints(n):
    return np.cos(np.arange(2*n + 4) * np.pi /(2*n + 3))

"Chebyshev collocation points in r"

col_r = SB_colpoints(PR)

colr = col_r[1:PR+2]

# colr1 = LR * colr/(np.sqrt(1-colr**2))                       # physical domain
rcol = LR * colr/(np.sqrt(1-colr**2))                       # physical domain 

# np.savetxt("rcolPR2L1.txt", rcol, fmt="%.16e", delimiter="\t")

# print("rcol", rcol)

" Base Matrix (Tchebyshev Polinomials) in r: "

SB_r = np.zeros([PR+1,PR+1])
dSB_r = np.zeros([PR+1,PR+1])
ddSB_r = np.zeros([PR+1,PR+1])

for i in range(PR+1):
  SB_r[i,] = SB(i,rcol)                                                 

for i in range(PR+1):
  dSB_r[i,] = dSB(i,rcol)

for i in range(PR+1):
  ddSB_r[i,] = ddSB(i,rcol)

" Big Matrix: "

Psi = np.kron(SB_r, P_x)

Psi_inv = np.linalg.inv(Psi)

drPsi = np.kron(dSB_r, P_x)

ddrPsi = np.kron(ddSB_r, P_x)

dxPsi = np.kron(SB_r, dP_x)

ddxPsi = np.kron(SB_r, ddP_x)

np.savetxt("PsiPR3L1PX3.txt", Psi, fmt="%.16e", delimiter="\t")
print("Psi saved ")

np.savetxt("Psi_invPR3L1PX3.txt", Psi_inv, fmt="%.16e", delimiter="\t")
print("Psi_inv saved ")

M = [[SB( i, rcol[k]) * P(2 * j, xcol[n])
          for i in range(PR + 1)
          for j in range(PX + 1)]
           for k in range(PR + 1)
           for n in range(PX + 1)]

M_inv = np.linalg.inv(M)
np.savetxt("MPR3L1PX3.txt", M, fmt="%.16e", delimiter="\t")
print("M saved ")

np.savetxt("M_invPR3L1PX3.txt", M_inv, fmt="%.16e", delimiter="\t")
print("M_inv saved ")