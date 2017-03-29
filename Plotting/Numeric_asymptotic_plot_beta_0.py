#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

# Asymptotic prediction for a range of beta values
beta_vals = [0]
#lines = ["k-","k--","k:","k-."]
lines = ["ks--", "k^:", "ko-."]
linecycler = cycle(lines)
for i in range(0, len(beta_vals) ) :
    # Load data
    data = np.loadtxt("./DATA/Parabolic_System/M_pi_K_beta_" + str(beta_vals[i]) + ".dat")
    K_asym = data[:,0]
    M_pi = -data[:,1]

    # Plot the asymptotic solution
    plt.plot( K_asym, M_pi, "k-")
    #plt.xlabel('K')
    #plt.ylabel('M / pi')

axes = plt.gca()
axes.set_xlim([0.0,1.2])
axes.set_ylim([-4.5,0.0])

# Numerical data for beta = 0
zeta_0_vals = [2, 4, 8]

# Load data
for i in range(0, len(zeta_0_vals) ) :
  #data2 = np.loadtxt("./DATA/K_" + str(K_num[i]) + "_beta_" + str(beta) + "_401x401_16_128" + "/A_file.dat")
  data2 = np.loadtxt("./DATA/K_Step_beta_0_zeta0_"+ str(zeta_0_vals[i]) + "/A_file.dat")
  zeta_0  = zeta_0_vals[i]
  K = data2[:,1]
  A = data2[:,2]
  # Plot the data (except the first data point)
  plt.plot( K[1:], A[1:] / zeta_0, next(linecycler), clip_on=False )

# Hide axis numbers
axes.xaxis.set_ticklabels([])
axes.yaxis.set_ticklabels([])

plt.savefig("Numeric_vs_asymptotic_beta_0.eps", format='eps', dpi=1000)
plt.show()
