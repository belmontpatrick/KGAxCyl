import numpy as np
import parameters as par
import basis as bs
import matplotlib.pyplot as plt
from matplotlib import cm
import collocation as col

A0 = par.A0
rcol = col.rcol
xcol = col.xcol
PX = par.PX
PR = par.PR
Basis = bs.Basis
P = bs.P
SB = bs.SB

"Importing Matrices:"
Psi = np.loadtxt("PsiPR3L1PX3.txt", delimiter="\t")

Psi_inv = np.loadtxt("PsiPR3L1PX3.txt", delimiter="\t")

M = np.loadtxt("MPR3L1PX3.txt", delimiter="\t")
M_inv = np.loadtxt("M_invPR3L1PX3.txt", delimiter="\t")

rtile = np.tile(rcol, (len(xcol)))

xrepeat = np.repeat(xcol, len(rcol))

Rcol, Xcol = np.meshgrid(rcol,xcol)

"Defining inital data"

def iniaxsph(A,r,x):
    return A * np.exp(-r**2)*(1 - x ** 2)

# inidata = iniaxsph(A0, rtile, xrepeat)
# phi0mesh = iniaxsph(A0, Rcol, Xcol)

# a0mesh = np.dot(M_inv, phi0mesh)

# a0sph = np.dot(Psi_inv, inidata)

phi_exact_col_arr = [iniaxsph(A0,col.rcol[j], col.xcol[k])
                          for j in range(PR + 1)
                          for k in range(PX + 1)]

a0_col_arr = np.dot(M_inv, phi_exact_col_arr)

#np.savetxt("a0sphPR3L1PX3.txt", a0sph, fmt="%.16e", delimiter="\t")

def phi_approx(a0,r,x):
    res = sum(a0[k] * Basis(k, r, x) for k in range((PR+1) * (PX+1)))
    return res

# phi0 = np.einsum('jl, l -> j', Psi_inv, phi_exact_col_arr)

Pcounts = (PR + 1) * (PX + 1)

da0sph = np.zeros(Pcounts)


####

# Pcounts = (PR + 1) * (PX + 1)

# ones = np.ones((PX+1))

# P_1 = np.zeros([PX+1,PX+1])

# for i in range(PX+1):
#   P_1[i,] = P(2*i,ones)

# zeros = np.zeros(rcol.shape)


# SB_0 = np.zeros([PR+1,PR+1])

# for i in range(PR+1):
#   SB_0[i,] = SB(i,zeros)

# print(SB_0.T[0])

# print(P_1.T[0])

# basis01 = np.kron(SB_0.T, P_1.T)

# b01 = basis01[0]

# # print(b01)

# phi001 = sum(a0sph[k] * b01[k] for k in range((PR+1) * (PX+1)))

# print(phi001)

rplot = np.linspace(0.00001, 1, 20)
xplot = np.linspace(-1, 1, 20)

Rplot, Xplot = np.meshgrid(rplot, xplot)

#phi0_approx = phi_approx(Rplot, Xplot, a0sph)
phi0_approx = phi_approx(a0_col_arr, Rplot, Xplot)

phiinidata = iniaxsph(A0, Rplot, np.flip(Xplot))

ax = plt.axes(projection = '3d')
# ax.plot_surface(Rplot, Xplot, phiinidata, cmap=cm.viridis)
ax.plot_surface(Rplot, Xplot, phi0_approx, cmap=cm.viridis)
plt.show()

# fig1 = plt.figure()
# ax1 = plt.axes(projection="3d")
# R, Z = np.meshgrid(np.linspace(0.0001, 10, 10), np.linspace(0.0001, 10, 10))
# Y1 = iniaxsph(A0,Rplot, Xplot)
# ax1.plot_surface(R, Z, Y1, cmap='plasma')
# ax1.set_xlabel('r')
# ax1.set_ylabel('z')
# ax1.set_zlabel('alpha0_exact')

# fig2 = plt.figure()
# ax2 = plt.axes(projection="3d")
# R, Z = np.meshgrid(np.linspace(0.0001, 10,0), np.linspace(0.0001, 10, 10))
# Y2 = phi_approx(R, Z, phi0)
# ax2.plot_surface(R, Z, Y2, cmap='plasma')
# ax2.set_xlabel('r')
# ax2.set_ylabel('z')
# ax2.set_zlabel('alpha0_approx')