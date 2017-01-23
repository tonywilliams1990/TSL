#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

beta = 0.1

# Load data
data = np.loadtxt("./DATA/Falkner_Skan_with_blowing_beta_" + str(beta) + "/eta_half_KB.dat")
KB  = data[:,0]
eta_half   = data[:,1]

# Plot the Falkner-Skan data
plt.plot(KB, eta_half, color='k')

plt.xlabel('KB')
axes = plt.gca()
axes.set_xlim([0,3])
axes.set_ylim([0,30])

# Extract the numerical data

K_num = [0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5]
eta_half_zeta_0_5 = []
eta_half_zeta_0_10 = []
eta_half_zeta_0_20 = []
# Load data
for i in range(0, len(K_num) ) :
  data2 = np.loadtxt("./DATA/K_" + str(K_num[i]) + "_beta_" + str(beta) + "_401x401_16_128" + "/A_file.dat")
  zeta_0  = data2[:,0]
  eta_half = data2[:,3]
  #print K[i]
  #print eta_half[19]
  #print zeta_0[ 19 ]
  eta_half_zeta_0_5.append( eta_half[4] )
  eta_half_zeta_0_10.append( eta_half[9] )
  eta_half_zeta_0_20.append( eta_half[19] )

# Plot the data
plt.plot( K_num, eta_half_zeta_0_5, 'ks--', clip_on=False )
plt.plot( K_num, eta_half_zeta_0_10, 'k^:', clip_on=False )
plt.plot( K_num, eta_half_zeta_0_20, 'ko-.', clip_on=False )

plt.show()
