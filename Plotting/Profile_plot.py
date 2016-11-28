#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt("./DATA/Parabolic_System/Base_soln.dat")
eta  = data[:,0]
UB   = data[:,1]
UBd  = data[:,2]
PhiB = data[:,3]
ThetaB = data[:,4]
ThetaBd = data[:,5]
PsiB = data[:,6]

#data_3D  = np.loadtxt("./DATA/Solution_3D.dat")
#beta_3D  = data_3D[:,0]
#shear_3D = data_3D[:,1]

# Plot the data
plt.plot(eta, UB, color='k', label="UB")
plt.plot(eta, UBd, linestyle='dashed', color='k', label='UBd' )
plt.plot(eta, PsiB, color='g', label='PsiB' )
plt.plot(eta, ThetaB, color='r', label='ThetaB' )


#plt.plot(beta_3D, shear_3D, linestyle='dashed' )
plt.xlabel('eta')
axes = plt.gca()
axes.set_xlim([0,30])
axes.set_ylim([0,1.5])
axes.legend()

plt.show()
