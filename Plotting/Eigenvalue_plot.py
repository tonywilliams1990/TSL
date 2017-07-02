#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

beta = 0

K_max = 0.8
K_min = -1.0
K_vals = np.linspace(K_min,K_max,19)
K_vals[10] = 0.0
#print K_vals

for i in range( 0, len( K_vals ) ):
  # Load data
  data = np.loadtxt("./DATA/Eigenvalue_problem_beta_" + str(beta) + "/Eigenvalues_K_" + str(K_vals[i]) + ".dat")
  #print K_vals[i]
  eigenvals = data
  for j in range( 0, eigenvals.size ):
    #print eigenvals[ j ]
    if eigenvals.size > 1:
      plt.plot( K_vals[i], eigenvals[j], 'ko', clip_on=False )
    else:
      plt.plot( K_vals[i], eigenvals, 'ko', clip_on=False )

axes = plt.gca()
axes.set_xlim([K_min,K_max])
axes.set_ylim([-1,3])

# Hide axis numbers
axes.xaxis.set_ticklabels([])
axes.yaxis.set_ticklabels([])

plt.savefig("Real_eigenvalue_beta_" + str(beta) + ".eps", format='eps', dpi=1000)

plt.show()
