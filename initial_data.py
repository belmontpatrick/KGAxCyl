import numpy as np
import basis as bs
import collocations as col
import parameters as par
import time
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

np.set_printoptions(precision=16)

t0 = time.process_time()

def alpha_m1_exact(r, z, A0, sigma_r, sigma_z):
    res = A0 * (np.exp( - (r**2 / sigma_r**2) - (z**2 / sigma_z**2) ) )
    return res

def alpha_m1_approx(r, z, a0, pr, pz):
    res = sum(a0[k] * bs.SBal(k, r, z) for k in range((pr+1) * (pz+1)))
    return res

PR_TEST = par.PR_TEST
PZ_TEST = par.PZ_TEST
LR_TEST = par.LR_TEST
LZ_TEST = par.LZ_TEST
A0_TEST = par.A0_TEST
R0_TEST = par.R0_TEST
SIGMA_R = par.SIGMA_R
SIGMA_Z = par.SIGMA_Z


alpha_m1_exact_col_arr = [alpha_m1_exact(col.col_r_shift[j], col.col_z_shift[k], A0_TEST, SIGMA_R, SIGMA_Z)
                         for j in range(PR_TEST + 1) for k in range(PZ_TEST + 1)]

alpha0 = np.einsum('jl,l->j', col.SBal_inv, alpha_m1_exact_col_arr)
dalpha0 = np.zeros((PR_TEST + 1)*(PZ_TEST + 1))

TRUNC = (PR_TEST + 1 )*(PZ_TEST +1)

fig1 = plt.figure()
ax1 = plt.axes(projection="3d")
r = np.linspace(0, 10, 100)
z = np.linspace(-5, 5, 100)
R, Z = np.meshgrid(r, z)
Y1 = alpha_m1_exact(R, Z, A0_TEST, SIGMA_R, SIGMA_Z)
ax1.plot_surface(R, Z, Y1, cmap = 'plasma')
ax1.set_xlabel('r')
ax1.set_ylabel('z')
ax1.set_zlabel('alpha0-1_exact')
# plt.savefig('alpha0_exact_cyl_trunc=%i.pdf' %(PR_TEST + 1)*(PZ_TEST + 1))
plt.savefig('alpha0-1_exact_cyl_trunc=%i.pdf' %TRUNC)

fig2 = plt.figure()
ax2 = plt.axes(projection="3d")
r = np.linspace(0, 10, 100)
z = np.linspace(-5, 5, 100)
R, Z = np.meshgrid(r, z)
Y2 = alpha_m1_approx(R, Z, alpha0, PR_TEST, PZ_TEST)
ax2.plot_surface(R, Z, Y2, cmap = 'plasma')
ax2.set_xlabel('r')
ax2.set_ylabel('z')
ax2.set_zlabel('alpha0-1_approx')
# plt.savefig('alpha0_approx_cyl_trunc=%i.pdf' %(PR_TEST + 1)*(PZ_TEST + 1))
plt.savefig('alpha0-1_approx_cyl_trunc=%i.pdf' %TRUNC)

fig3 = plt.figure()
ax3 = plt.axes(projection="3d")
error = alpha_m1_exact(R, Z, A0_TEST, SIGMA_R, SIGMA_Z) - alpha_m1_approx(R, Z, alpha0, PR_TEST, PZ_TEST)
ax3.plot_surface(R, Z, error, cmap = 'plasma')
ax3.set_xlabel('r')
ax3.set_ylabel('z')
ax3.set_zlabel('error')
# plt.savefig('alpha0_approx_cyl_trunc=%i.pdf' %(PR_TEST + 1)*(PZ_TEST + 1))
plt.savefig('error_trunc=%i.pdf' %TRUNC)

plt.show()

t1 = time.process_time()
print('Running time (initial data):', t1-t0, 's')