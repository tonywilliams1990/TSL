#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

K = 2
beta = 0.1

# Load data
data = np.loadtxt("./DATA/Falkner_Skan_with_blowing_beta_0.1/Base_soln_KB_2.dat")
eta  = data[:,0]
UB   = data[:,1]
UBd  = data[:,2]
PhiB = data[:,3]
ThetaB = data[:,4]
ThetaBd = data[:,5]
PsiB = data[:,6]


# Plot the data
plt.plot(eta, UB, color='k', label="U Falkner beta=0.1, K=2")
#plt.plot(eta, UBd, linestyle='dashed', color='k', label='UBd' )
#plt.plot(eta, PsiB, color='g', label='PsiB' )
#plt.plot(eta, ThetaB, color='r', label='ThetaB' )

# Load numerical data
zeta0s = [20, 40]
lines = ["-","--","-.",":"]
linecycler = cycle(lines)

for i in range(0, len(zeta0s) ) :
  data2 = np.loadtxt("./DATA/K_" + str(K) + "_beta_" + str(beta) + "_401x401_16_128" + "/Qout_" + str(zeta0s[i]) + ".dat")
  eta = data2[:,1]  
  U  = data2[:,8]
  eta_hzeta_0 = eta[0:401]
  U_hzeta_0 = U[0:401]

  plt.plot(eta_hzeta_0, U_hzeta_0, next(linecycler), label="U_num zeta0=" + str(zeta0s[i]) )


# axes etc
plt.xlabel('eta')
axes = plt.gca()
axes.set_xlim([0,30])
axes.set_ylim([0,1.5])
axes.legend()

plt.savefig("profile_plot.eps", format='eps', dpi=1000)

plt.show()
