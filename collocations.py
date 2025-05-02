import time
import numpy as np
import basis as bs
import parameters as par
# import matplotlib.pyplot as plt

t0 = time.process_time()

np.set_printoptions(precision=16)


PR_TEST = par.PR_TEST
PZ_TEST = par.PZ_TEST
LR_TEST = par.LR_TEST
LZ_TEST = par.LZ_TEST


def colx(k):
    res = np.cos(np.pi * k / (2 * PR_TEST + 3))
    return res


def coly(k):
    res = np.cos(np.pi * k / (2 * PZ_TEST + 3))
    return res


def colr(k):
    res = LR_TEST * colx(k) / (np.sqrt(1 - (colx(k)**2)))
    return res


def colz(k):
    res = LZ_TEST * coly(k) / (np.sqrt(1 - (coly(k)**2)))
    return res


"""Collocation points shifted"""

PR_SHIFT = PR_TEST + 1
PZ_SHIFT = PZ_TEST + 1

col_r_shift = np.zeros(PR_SHIFT)

for j in range(0, PR_SHIFT):
    col_r_shift[j] = colr(j + 1)

# print(col_r_shift)

col_z_shift = np.zeros(PZ_SHIFT)

for j in range(0, PZ_SHIFT):
    col_z_shift[j] = colr(j + 1)

# print(col_z_shift)

# r_teste = np.zeros((PR_TEST + 1, PR_TEST + 1))

# for j in range(PR_TEST + 1):
#     for k in range(PR_TEST + 1):
#         r_teste[j][k] = col_r_shift[j] 

# r = [ r_teste[j][k] for j in range(PR_TEST + 1) for k in range(PZ_TEST + 1)]



'''Matrices for alpha basis at pr+1 and pz+1 collocation points'''

SBal = [[bs.SBr(2 * i, col_r_shift[k]) * bs.SBz(2 * j, col_z_shift[n])
            for i in range(PR_TEST + 1)
            for j in range(PZ_TEST + 1)]
           for k in range(PR_TEST + 1)
           for n in range(PZ_TEST + 1)]

drSBal = [[bs.dSBr(2 * i, col_r_shift[k]) * bs.SBz(2 * j, col_z_shift[n])
            for i in range(PR_TEST + 1)
            for j in range(PZ_TEST + 1)]
           for k in range(PR_TEST + 1)
           for n in range(PZ_TEST + 1)]

ddrSBal = [[bs.ddSBr(2 * i, col_r_shift[k]) * bs.SBz(2 * j, col_z_shift[n])
            for i in range(PR_TEST + 1)
            for j in range(PZ_TEST + 1)]
           for k in range(PR_TEST + 1)
           for n in range(PZ_TEST + 1)]

dzSBal = [[bs.SBr(2 * i, col_r_shift[k]) * bs.dSBz(2 * j, col_z_shift[n])
            for i in range(PR_TEST + 1)
            for j in range(PZ_TEST + 1)]
           for k in range(PR_TEST + 1)
           for n in range(PZ_TEST + 1)]

ddzSBal = [[bs.SBr(2 * i, col_r_shift[k]) * bs.ddSBz(2 * j, col_z_shift[n])
            for i in range(PR_TEST + 1)
            for j in range(PZ_TEST + 1)]
           for k in range(PR_TEST + 1)
           for n in range(PZ_TEST + 1)]

SBal_inv = np.linalg.inv(SBal)



"Kron way"

# " Base Matrix (Tchebyshev Polinomials) in z: "

# SB_z = np.zeros([PZ+1,PZ+1])
# zSB = np.zeros([PZ+1,PZ+1])
# zzSB = np.zeros([PZ+1,PZ+1])

# for i in range(PZ+1):
#   SB_z[i,] = bs.SBz(i,col_z_shift)                                                 

# for i in range(PZ+1):
#   zSB[i,] = bs.dSBz(i,col_z_shift)

# for i in range(PZ+1):
#   zzSB[i,] = bs.ddSBz(i,col_z_shift)

# # np.savetxt("SB_zPZ30L5.txt", SB_z, fmt="%.16e", delimiter="\t")

# " Base Matrix (Tchebyshev Polinomials) in r: "

# SB_r = np.zeros([PR+1,PR+1])
# rSB = np.zeros([PR+1,PR+1])
# rrSB = np.zeros([PR+1,PR+1])

# for i in range(PR+1):
#   SB_r[i,] = bs.SBr(i,col_r_shift)                                                 

# for i in range(PR+1):
#   rSB[i,] = bs.dSBr(i,col_r_shift)

# for i in range(PR+1):
#   rrSB[i,] = bs.ddSBr(i,col_r_shift)                     


# SB = np.kron(SB_r.T,SB_z.T)

# SB_inv = np.linalg.inv(SB)

# drSB = np.kron(rSB.T,SB_z.T)

# ddrSB = np.kron(rrSB.T,SB_z.T)

# dzSB = np.kron(SB_r.T,zSB.T)

# ddzSB = np.kron(SB_r.T,zzSB.T)

# print("SB_inv")
# print(SB_inv)

t1 = time.process_time()
print('Running time (collocationpoints):', t1 - t0, 's')