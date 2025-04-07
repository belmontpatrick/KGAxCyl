import numpy as np

#defining grid in time

t0 = 0.0
tf = 1.0
delta = t0 - tf
h = 0.1
N = int((tf - t0) / h)
#t = np.linspace(t0,tf,N)