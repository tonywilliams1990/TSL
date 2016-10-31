#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("./DATA/Solution_2D.dat")
beta  = data[:,0]
shear = data[:,1]

data_3D  = np.loadtxt("./DATA/Solution_3D.dat")
beta_3D  = data_3D[:,0]
shear_3D = data_3D[:,1]

# Plot the data
plt.plot(beta, shear, color='k')
plt.plot(beta_3D, shear_3D, linestyle='dashed' )
plt.xlabel('beta')
plt.ylabel('U_B\'(eta=0)')
axes = plt.gca()
axes.set_xlim([-0.5,1.0])
axes.set_ylim([-0.5,1.25])

plt.show()
