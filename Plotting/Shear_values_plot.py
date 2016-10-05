#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("./DATA/Shear_values.dat")
beta  = data[:,0]
shear = data[:,1]

# Plot the data
plt.plot(beta, shear, color='k')
plt.xlabel('beta')
plt.ylabel('U_B\'(eta=0)')
axes = plt.gca()
axes.set_xlim([-0.45,1.0])
axes.set_ylim([-0.2,1.2])

plt.show()
