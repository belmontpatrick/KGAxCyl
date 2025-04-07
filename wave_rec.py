import basis as bs
import numpy as np
import parameters as pr

PR = pr.PR

line_0 = 0

p_test = 10

acyl_datatest = np.loadtxt("rk4_acyl.txt")[0,0]
print(acyl_datatest.shape)

# a00 = np.loadtxt("rk4_acyl.txt")[line_0, 1:PR + 2]

# #a0 = np.loadtxt("runge_kutta_BSSN.txt")[line_0, 1:p_test+2]

# print(a00)

