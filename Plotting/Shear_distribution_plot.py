#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

K_val = [0.5, 1, 2]
beta = 0.1
M = 400
N = 400
hzeta_inf = 16
eta_inf = 128

for i in range(0, len(K_val) ) :

  ########################### Plot 1 #############################
  zeta0 = 20
  # Load data
  data = np.loadtxt("./DATA/K_" + str(K_val[i]) + "_beta_" + str(beta) + "_" + str(N+1) + "x" + str(M+1) + "_" + str(hzeta_inf) + "_" + str(eta_inf) + "/Wall_shear_zeta0_" + str(zeta0) + ".dat")
  hzeta  = data[:,0]
  shear = data[:,1]

  # Plot the data
  plt.plot(hzeta, shear, color='k', linestyle='dashed')
  #plt.xlabel('hzeta')
  #plt.ylabel('U_B\'(eta=0)')

  ############################ Plot 2 ############################
  zeta0 = 40
  # Load data
  data2 = np.loadtxt("./DATA/K_" + str(K_val[i]) + "_beta_" + str(beta) + "_" + str(N+1) + "x" + str(M+1) + "_" + str(hzeta_inf) + "_" + str(eta_inf) + "/Wall_shear_zeta0_" + str(zeta0) + ".dat")
  hzeta2  = data2[:,0]
  shear2 = data2[:,1]

  # Plot the data
  plt.plot(hzeta2, shear2, color='k', linestyle=':')

  ######################## Parabolic plot ########################
  base = "2D"
  # Load data
  data3 = np.loadtxt("./DATA/Parabolic_System/Wall_shear_K_" + str(K_val[i]) + "_beta_" + str(beta) + "_" + base + ".dat")
  hzeta3  = data3[:,0]
  shear3 = data3[:,1]

  # Plot the data (thick line)
  plt.plot(hzeta3, shear3, color='k', linewidth=1.0)

axes = plt.gca()
axes.set_xlim([0.0,2.0])
axes.set_ylim([0.0,0.6])

# Hide axis numbers
axes.xaxis.set_ticklabels([])
axes.yaxis.set_ticklabels([])

plt.savefig("Shear_distribution_beta_" + str(beta) + ".eps", format='eps', dpi=1000)
plt.show()
