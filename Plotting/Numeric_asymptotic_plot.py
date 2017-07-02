#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

# Asymptotic prediction for a range of beta values
beta_vals = [0.05, 0.1, 0.4, 0.8]
#lines = ["k-","k--","k:","k-."]
lines = ["k-"]
linecycler = cycle(lines)
for i in range(0, len(beta_vals) ) :
    # Load data
    data = np.loadtxt("./DATA/Parabolic_System/M_pi_K_beta_" + str(beta_vals[i]) + ".dat")
    K_asym = data[:,0]
    M_pi = -data[:,1]

    # Plot the asymptotic solution
    plt.plot( K_asym, M_pi, next(linecycler))
    #plt.xlabel('K')
    #plt.ylabel('M / pi')

axes = plt.gca()
axes.set_xlim([0.0,3.0])
axes.set_ylim([-10.0,0.0])

# Numerical data for beta = 0.1
beta = 0.1
K_num = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3]
A_zeta_0_2 = []
A_zeta_0_4 = []
A_zeta_0_8 = []
# Load data
for i in range(0, len(K_num) ) :
  data2 = np.loadtxt("./DATA/K_" + str(K_num[i]) + "_beta_" + str(beta) + "_401x401_16_128" + "/A_file.dat")
  zeta_0  = data2[:,0]
  A = data2[:,1]
  #print K[i]
  #print A[1]
  #print zeta_0[ 1 ]
  A_zeta_0_2.append( A[1] / zeta_0[ 1 ] )
  A_zeta_0_4.append( A[3] / zeta_0[ 3 ] )
  A_zeta_0_8.append( A[7] / zeta_0[ 7 ] )

#print A_zeta_0_2
#print A_zeta_0_4
#print A_zeta_0_8

# Plot the data
plt.plot( K_num, A_zeta_0_2, 'ks--', clip_on=False )
plt.plot( K_num, A_zeta_0_4, 'k^:', clip_on=False )
plt.plot( K_num, A_zeta_0_8, 'ko-.', clip_on=False )

# Hide axis numbers
axes.xaxis.set_ticklabels([])
axes.yaxis.set_ticklabels([])

plt.savefig("Numeric_vs_asymptotic_beta_" + str(beta) + ".eps", format='eps', dpi=1000)
plt.show()
