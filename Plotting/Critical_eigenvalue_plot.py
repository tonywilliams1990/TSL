#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

K = [ -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6 ]
beta = [ -0.2917, -0.0187, 0.1944, 0.3473, 0.4419, 0.4862, 0.4990, 0.5009, 0.5011, 0.5012, 0.5013, 0.5014, 0.5015, 0.5016, 0.5017]

N_interp = 1000


# Interpolate the data
#eta_new = np.linspace( eta[0], eta[-1], N_interp, endpoint=True)
#u_new   = interp1d( eta, u, kind='cubic' )
#psi_new = interp1d( eta, psi, kind='cubic' )
K_new = np.linspace( K[0], K[-1], N_interp, endpoint=True)
beta_new = interp1d( K, beta, kind='cubic' )

# Plot data
plt.plot( K, beta, 'ko')
#plt.plot( eta_new, u_new(eta_new), 'k', label = 'u' )
plt.plot( K_new, beta_new(K_new), 'k' )

# Plot Rich's data
data = np.loadtxt("./DATA/nc1.dat")
beta_rich  = data[:,0]
K_rich = data[:,1]
plt.plot( -K_rich, beta_rich, 'r' )



axes = plt.gca()
axes.set_xlim([-1,6])
axes.set_ylim([-0.3,0.6])
plt.xlabel('K')
plt.ylabel('beta')
plt.title('lambda_r = 0')

# Hide axis numbers
#axes.xaxis.set_ticklabels([])
#axes.yaxis.set_ticklabels([])
#plt.legend()

plt.savefig("Critcal_eigenvalue.eps", format='eps', dpi=1000)

plt.show()
