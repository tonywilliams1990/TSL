#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

beta = 0.1

K = 2.0


# Load data
data = np.loadtxt("./DATA/Eigenvalue_problem_beta_" + str(beta) + "/Eigenvectors_K_" + str(K) + ".dat")

eta = data[:,0]
u   = data[:,1]

u_max = np.max( u )
u_min = np.min( u )
min_index = np.where( u == u_min )
max_index = np.where( u == u_max )

print eta[min_index]
print eta[max_index]

# Plot data
plt.plot( eta, u, 'k' )

axes = plt.gca()
axes.set_xlim([0,30])
axes.set_ylim([u_min,u_max])

# Hide axis numbers
#axes.xaxis.set_ticklabels([])
#axes.yaxis.set_ticklabels([])

plt.savefig("eigenvectors_beta_" + str(beta) + ".eps", format='eps', dpi=1000)

plt.show()
