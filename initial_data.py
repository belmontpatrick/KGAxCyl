import numpy as np
import collocations as col
import parameters as par

A0 = par.A0
PR = par.PR
PZ = par.PZ
rcol = col.rcol
zcol = col.zcol
SB = col.SB
SB_inv = col.SB_inv

def cylindricalgaussian(A,r,z):

    return A * np.exp(-r**2 - z**2)

##tiling rcol

rtile = np.tile(rcol, (len(zcol)))

#repeating zcol
zrep = np.repeat(zcol, len(rcol))

#puting initial data

g0cyl = cylindricalgaussian(A0,rtile,zrep)

#obtaning coeficient in t = 0

a0cyl = np.dot(SB_inv,g0cyl)

#defining da0cyl

da0cyl = np.zeros(((PR + 1) * (PZ + 1),1))

