import numpy as np
from initial_data import alpha0, dalpha0
from timegrid import h, It, t
from dynsys import dda, phi
import time


t0 = time.process_time()

np.set_printoptions(precision=16)


with open('resultados_a_def_tst.txt', 'w') as a0_file, \
     open('resultados_da_def_tst.txt', 'w') as da_file:

    for i in range(It):  # Runge Kutta 4th order

        phi_ = phi(alpha0)
        dda_ = dda(alpha0)
        L1 = h*(dalpha0)
        K1 = h*(dda_)

        phi_ = phi(alpha0 + L1/2)
        dda_ = dda(alpha0 + L1/2)
        L2 = h*(dalpha0 + K1/2)
        K2 = h*(dda_)

        phi_ = phi(alpha0 + L2/2)
        dda_ = dda(alpha0 + L2/2)
        L3 = h*(dalpha0 + K2/2)
        K3 = h*(dda_)

        phi_ = phi(alpha0 + L3)
        dda_ = dda(alpha0 + L3)
        L4 = h*(dalpha0 + K3)
        K4 = h*(dda_)

        dalpha0 = dalpha0 + 1/6 * (K1 + 2*K2 + 2*K3 + K4)
        alpha0 = alpha0 + 1/6 * (L1 + 2*L2 + 2*L3 + L4)
        
        a0_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in alpha0.flatten()])
        da_line = f"{t[i]:.6f} " + " ".join([f"{val:.15f}" for val in dalpha0.flatten()])
        
        a0_file.write(a0_line + "\n")
        da_file.write(da_line + "\n")


t1 = time.process_time()
print('Running time (rk4):', t1 - t0, 's')