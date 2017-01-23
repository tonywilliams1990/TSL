#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("./DATA/Solution_test.dat")
x  = data[:,0]
y   = data[:,1]
yd  = data[:,2]

#data_3D  = np.loadtxt("./DATA/Solution_3D.dat")
#beta_3D  = data_3D[:,0]
#shear_3D = data_3D[:,1]

# Plot the data
plt.plot(x, y, color='k', label="y(x)")
plt.plot(x, 0.5*yd*yd - np.cos(y), color='r', label="0.5*y'(x)^2 - cos(y)")
plt.plot(x, 0.25*yd*yd, color='g', label="c_0*|y'(x)|^2")
plt.plot(x, 0.5*yd*yd + 1, color='b', label="c_1*|y'(x)|^2 + c_2")

#plt.plot(beta_3D, shear_3D, linestyle='dashed' )
plt.xlabel('x')
axes = plt.gca()
axes.set_xlim([0,1])
axes.set_ylim([0,10])
axes.legend()

plt.show()
