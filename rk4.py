import numpy as np
import basis as bs
import collocations as col
import initial_data as ini
import timegrid as gr
import dynsys as ds
import plot_matrices as pm
import parameters as par

np.set_printoptions(precision=16)

PR = par.PR
PZ = par.PZ
M = pm.M
N = gr.N
SB_plot = pm.SB_plot
a0cyl = ini.a0cyl
da0cyl = ini.da0cyl
h = gr.h
#t = gr.t
N = gr.N
phi = ds.phi
dda = ds.dda

#Runge-Kutta 4th order


file_acyl = open("rk4_acyl.txt", "w")
# file_dacyl = open("rk4_dacyl.txt", "w")

t = 0.0           # Initial time

for i in range(N):  # Runge Kutta 4th order

  phi_ = phi(a0cyl)
  dda_ = dda(a0cyl)
  L1 = h*da0cyl
  K1 = h*(dda_)

  phi_ = phi(a0cyl + L1/2)
  dda_ = dda(a0cyl + L1/2)
  L2 = h*(da0cyl + K1/2)
  K2 = h*(dda_)

  phi_ = phi(a0cyl + L2/2)
  dda_ = dda(a0cyl + L2/2)
  L3 = h*(da0cyl + K2/2)
  K3 = h*(dda_)

  phi_ = phi(a0cyl + L3)
  dda_ = dda(a0cyl + L3)
  L4 = h*(da0cyl + K3)
  K4 = h*(dda_)

  da0cyl = da0cyl + 1/6 * (K1 + 2*K2 + 2*K3 + K4)
  a0cyl = a0cyl + 1/6 * (L1 + 2*L2 + 2*L3 + L4)

  t += h
  # Save the results to the files
  file_acyl.write(str(t))

  for j in range(len(a0cyl)):
              file_acyl.write(str(a0cyl[j]))
              file_acyl.write(' ')

  file_acyl.write("\n")

# Close files
file_acyl.close()
# file_dacyl.close()

