#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

beta_max = 0.7
beta_vals = np.linspace(0.1, beta_max, 7)

K = 1.5

for i in range( 0, len( beta_vals ) ):
  # Load data
  if beta_vals[i] == 0.0:
    data = np.loadtxt("./DATA/Eigenvalue_problem_beta_" + str(0) + "/Eigenvalues_K_" + str(K) + ".dat")
  else:
    data = np.loadtxt("./DATA/Eigenvalue_problem_beta_" + str(beta_vals[i]) + "/Eigenvalues_K_" + str(K) + ".dat")
  #print K_vals[i]
  eigenvals = data
  for j in range( 0, eigenvals.size ):
    #print eigenvals[ j ]
    if eigenvals.size > 1:
      plt.plot( beta_vals[i], eigenvals[j], 'ko', clip_on=False )
    else:
      plt.plot( beta_vals[i], eigenvals, 'ko', clip_on=False )

axes = plt.gca()
axes.set_xlim([0,beta_max])
axes.set_ylim([-1,3])

# Hide axis numbers
#axes.xaxis.set_ticklabels([])
#axes.yaxis.set_ticklabels([])

plt.savefig("Real_eigenvalue_K_" + str(K) + ".eps", format='eps', dpi=1000)

plt.show()
