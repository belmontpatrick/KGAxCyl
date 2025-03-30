import numpy as np
import basis as bs
import repelem as rep
import matplotlib.pyplot as plt
import parameters as par
import collocations as col

SBz = bs.SBz
SBr = bs.SBr
dSBz = bs.dSBz
dSBr = bs.dSBr
ddSBr = bs.ddSBr
ddSBz = bs.ddSBz
repelem = rep.repelem
PR = par.PR
PZ = par.PZ
SB = col.SB

M = 100

rplot = np.linspace(0.000001,10,M)
zplot = np.linspace(0.000001,10,M)

# SBphi_plot =  [[SBr(2 * i, rplot[k]) * SBz(2 * j, zplot[n])
#             for i in range(PR + 1)
#             for j in range(PZ + 1)]
#            for k in range(PR + 1)
#            for n in range(PZ + 1)]

Rplot, Zplot = np.meshgrid(rplot, zplot, indexing = 'ij') 

SBz_plot = np.zeros([PZ+1,M])   ###PZ + 1, M
SBr_plot = np.zeros([PR+1,M])   ###PR + 1, M
dSBz_plot = np.zeros([PZ+1,M])
dSBr_plot = np.zeros([PR+1,M])
ddSBz_plot = np.zeros([PZ+1,M])
ddSBr_plot = np.zeros([PR+1,M])

for i in range(PR+1):
  SBr_plot[i,] = SBr(i,rplot) 

for i in range(PZ+1):
  SBz_plot[i,] = SBr(i,zplot)

for i in range(PR+1):
  dSBr_plot[i,] = dSBr(i,rplot) 

for i in range(PR+1):
  ddSBr_plot[i,] = SBr(i,zplot) 

for i in range(PZ+1):
  dSBz_plot[i,] = SBr(i,zplot)

for i in range(PZ+1):
  ddSBz_plot[i,] = SBr(i,zplot)

# psirplot = SBr_plot[]

SB_plot = np.tile(SBz_plot,(PR+1,M)) * repelem(SBr_plot,(PZ+1,M))