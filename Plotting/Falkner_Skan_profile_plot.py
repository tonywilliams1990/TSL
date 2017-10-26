#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

#K = 2
#beta = 0.1
K_vals = [0,0.2,0.4,0.6,0.8]

for K in range(0, len(K_vals) ) :

    # Load data
    data = np.loadtxt("./DATA/Falkner_Skan_with_blowing_beta_0/Base_soln_KB_" + str(K_vals[K]) + ".dat")
    eta  = data[:,0]
    UB   = data[:,1]
    UBd  = data[:,2]
    PhiB = data[:,3]
    ThetaB = data[:,4]
    ThetaBd = data[:,5]
    PsiB = data[:,6]

    # Plot the data
    #plt.figure()
    plt.plot(UB + K + 0.2*K, eta, color='k')
    #plt.plot(eta, UBd, linestyle='dashed', color='k', label='UBd' )
    #plt.plot(eta, PsiB, color='g', label='PsiB' )
    #plt.plot(eta, ThetaB, color='r', label='ThetaB' )

    # axes etc
    #plt.xlabel('eta')
    axes = plt.gca()
    axes.set_xlim([0,6])
    axes.set_ylim([0,12])
    #axes.legend()

    # Hide axis numbers
    axes.axis('off')
    #axes.xaxis.set_ticklabels([])
    #axes.yaxis.set_ticklabels([])

plt.savefig("profile_plots.eps", format='eps', dpi=1000)

plt.show()
