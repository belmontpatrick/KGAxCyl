import numpy as np
import parameters as par
import basis as bs
import repelem as rep

np.set_printoptions(precision=16)

LR = par.LR
PR = par.PR
PZ = par.PZ
LZ = par.LZ
LR = par.LR
SBr = bs.SBr
dSBr = bs.dSBr
ddSBr = bs.ddSBr
SBz = bs.SBz
dSBz = bs.dSBz
ddSBz = bs.ddSBz
repelem = rep.repelem

"Chebyshev collocation points"

def SB_colpoints(n):
    return np.cos(np.arange(2*n + 4)*np.pi /(2*n + 3))

"Chebyshev collocation points in z"
col_z = SB_colpoints(PZ)

colz = col_z[1:PZ+2]

z1 = LZ * colz/(np.sqrt(1-colz**2))                       # physical domain
zcol = np.flip(z1)

#np.savetxt("zcolPZ30L5.txt", zcol, fmt="%.16e", delimiter="\t")

"Chebyshev collocation points in r"

col_r = SB_colpoints(PR)

colr = col_r[1:PR+2]

colr1 = LR * colr/(np.sqrt(1-colr**2))                       # physical domain 
rcol = np.flip(colr1)

# np.savetxt("rcolPR60L5.txt", rcol, fmt="%.16e", delimiter="\t")


" Base Matrix (Tchebyshev Polinomials) in z: "

SB_z = np.zeros([PZ+1,PZ+1])
zSB = np.zeros([PZ+1,PZ+1])
zzSB = np.zeros([PZ+1,PZ+1])

for i in range(PZ+1):
  SB_z[i,] = SBz(i,zcol)                                                 

for i in range(PZ+1):
  zSB[i,] = dSBz(i,zcol)

for i in range(PZ+1):
  zzSB[i,] = ddSBz(i,zcol)

# np.savetxt("SB_zPZ30L5.txt", SB_z, fmt="%.16e", delimiter="\t")

" Base Matrix (Tchebyshev Polinomials) in r: "

SB_r = np.zeros([PR+1,PR+1])
rSB = np.zeros([PR+1,PR+1])
rrSB = np.zeros([PR+1,PR+1])

for i in range(PR+1):
  SB_r[i,] = SBr(i,rcol)                                                 

for i in range(PR+1):
  rSB[i,] = dSBr(i,rcol)

for i in range(PR+1):
  rrSB[i,] = ddSBr(i,rcol)                     

#np.savetxt("SB_rPZ60L5.txt", SB_r, fmt="%.16e", delimiter="\t")

" Cylindrical Basis Matrix :"

SB = np.kron(SB_r,SB_z)

SB_inv = np.linalg.inv(SB)

drSB = np.kron(rSB,SB_z)

ddrSB = np.kron(rrSB,SB_z)

dzSB = np.kron(SB_r,zSB)

ddzSB = np.kron(SB_r,zzSB)
