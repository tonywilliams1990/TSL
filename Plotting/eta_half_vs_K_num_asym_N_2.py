#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

beta = 0
zeta0 = 20

# Load data
data = np.loadtxt("./DATA/K_Step_beta_" + str(beta) + "_zeta0_" + str(zeta0) + "_N_2/A_file.dat")
K  = data[:,1]
eta_half   = data[:,4]
K_max = np.max(K)

# Plot the numerical data
plt.plot(K, eta_half, color='k')

plt.xlabel('K')
axes = plt.gca()
axes.set_xlim([0,4])
axes.set_ylim([0,20])

# Plot the asymptotic predictions
F_minus_inf = -1.238
K_crit_1 = -np.pi * F_minus_inf / 2

K_1 = np.linspace( K_crit_1, K_max)
r_s_1 = ( - 2 * K_1 * zeta0 ) / ( F_minus_inf * np.pi - 2 * K_1 )
plt.plot(K_1, r_s_1, color='r')


#plt.savefig("plot.eps", format='eps', dpi=1000)

plt.show()
