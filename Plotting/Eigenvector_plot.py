#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

beta = 0.1
K = 2.0
N_interp = 1000


# Load data
data = np.loadtxt("./DATA/Eigenvalue_problem_beta_" + str(beta) + "/Eigenvectors_K_" + str(K) + ".dat")

eta    = data[:,0]
u      = data[:,1]
ud     = data[:,2]
phi    = data[:,3]
theta  = data[:,4]
thetad = data[:,5]
psi    = data[:,6] 

u_max = np.max( u )
u_min = np.min( u )
min_index = np.where( u == u_min )
max_index = np.where( u == u_max )

print eta[min_index]
print eta[max_index]

# Interpolate the data
eta_new = np.linspace( eta[0], eta[-1], N_interp, endpoint=True)
u_new   = interp1d( eta, u, kind='cubic' )
psi_new = interp1d( eta, psi, kind='cubic' )

# Plot data
#plt.plot( eta_new, u_new(eta_new), 'k', label = 'u' )
plt.plot( eta_new, psi_new(eta_new), 'r-', label = 'psi' )
plt.plot( eta_new, (1-beta)*u_new(eta_new), 'k', label = '(1-beta)*u' )

axes = plt.gca()
axes.set_xlim([0,30])
#axes.set_ylim([u_min,u_max])
axes.set_ylim([-0.05,0.15])

# Hide axis numbers
#axes.xaxis.set_ticklabels([])
#axes.yaxis.set_ticklabels([])
plt.legend()

plt.savefig("eigenvectors_beta_" + str(beta) + "_K_" + str(K) + ".eps", format='eps', dpi=1000)

plt.show()
