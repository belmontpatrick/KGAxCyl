# import numpy as np

# PR = 20
# PZ = 10
# LR = 1
# LZ = 1
# A0 = 0.0002

import time
import numpy as np


np.set_printoptions(precision=16)

t0 = time.process_time()

PR_TEST = 10
PZ_TEST = 10
LR_TEST = 1.0
LZ_TEST = 1.0
A0_TEST = 0.002
R0_TEST = 5.0
SIGMA_R = 1.0
SIGMA_Z = 1.0
T = 1.0
N = 1000

t1 = time.process_time()
print('Running time (parameters):', t1 - t0, 's')