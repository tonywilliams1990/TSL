#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

beta = 0
K = 2
# change the levels in the VW plot to see more detail (line 85)

data = np.loadtxt("./DATA/Marching_K_"+ str(K) + "_beta_" + str(beta) + "_1000x201x201_100_30_30/A_file.dat")

zeta0 = data[:,0]
A = data[:,1]
x = data[:,2]
shear_at_origin = data[:,3]
streak_radius = data[:,4]
integral_U_squared = data[:,5]

# A plot
plt.plot(x, A, 'k-')
plt.savefig('A_vs_x_K_' + str(K) + '_beta_' + str(beta) + '.eps', format='eps', dpi=1000)

# Shear plot
plt.figure()
plt.plot(x, shear_at_origin, 'k-')
plt.savefig('shear_at_origin_vs_x_K_' + str(K) + '_beta_' + str(beta) + '.eps', format='eps', dpi=1000)

# integral U^2 plot
plt.figure()
plt.plot(x, integral_U_squared, 'k-')
plt.savefig('integral_U_squared_vs_x_K_' + str(K) + '_beta_' + str(beta) + '.eps', format='eps', dpi=1000)



plt.show()
