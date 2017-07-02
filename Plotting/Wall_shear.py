#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

K = 2.5
beta = 0
M = 400
N = 400
hzeta_inf = 16
eta_inf = 128

########################### Plot 1 #############################
zeta0 = 20
# Load data
data = np.loadtxt("./DATA/K_" + str(K) + "_beta_" + str(beta) + "_" + str(N+1) + "x" + str(M+1) + "_" + str(hzeta_inf) + "_" + str(eta_inf) + "/Wall_shear_zeta0_" + str(zeta0) + ".dat")
hzeta  = data[:,0]
shear = data[:,1]

# Plot the data
plt.plot(hzeta, shear, color='r', linestyle='dashed')
plt.xlabel('hzeta')
plt.ylabel('U_B\'(eta=0)')
axes = plt.gca()
axes.set_xlim([0.0,2.0])
axes.set_ylim([0.0,0.5])

############################ Plot 2 ############################
zeta0 = 40
# Load data
data2 = np.loadtxt("./DATA/K_" + str(K) + "_beta_" + str(beta) + "_" + str(N+1) + "x" + str(M+1) + "_" + str(hzeta_inf) + "_" + str(eta_inf) + "/Wall_shear_zeta0_" + str(zeta0) + ".dat")
hzeta2  = data2[:,0]
shear2 = data2[:,1]

# Plot the data
plt.plot(hzeta2, shear2, color='r')

######################## Parabolic plot ########################
base = "2D"
# Load data
data3 = np.loadtxt("./DATA/Parabolic_System/Wall_shear_K_" + str(K) + "_beta_" + str(beta) + "_" + base + ".dat")
hzeta3  = data3[:,0]
shear3 = data3[:,1]

# Plot the data
plt.plot(hzeta3, shear3, color='k')

plt.show()
