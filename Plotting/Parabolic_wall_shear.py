#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

K = 0.5
beta = 0.1
base = "3D"

# Load data
data = np.loadtxt("./DATA/Parabolic_System/Wall_shear_K_" + str(K) + "_beta_" + str(beta) + "_" + base + ".dat")
hzeta  = data[:,0]
shear = data[:,1]

# Plot the data
plt.plot(hzeta, shear, color='k')
plt.xlabel('hzeta')
plt.ylabel('U_B\'(eta=0)')
axes = plt.gca()
axes.set_xlim([0.0,2.0])
axes.set_ylim([0.0,1.0])

plt.show()
