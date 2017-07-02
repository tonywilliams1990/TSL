#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

#K = 2
beta_vals = [0.1,0.5, 0.8]
K = 9
plot_full_data = False

for beta in range(0, len(beta_vals) ) :

    # Load data
    data = np.loadtxt("./DATA/Falkner_Skan_with_blowing_beta_" + str(beta_vals[beta]) + "/Base_soln_KB_" + str(K) + ".dat")
    eta  = data[:,0]
    UB   = data[:,1]
    UBd  = data[:,2]
    PhiB = data[:,3]
    ThetaB = data[:,4]
    ThetaBd = data[:,5]
    PsiB = data[:,6]

    # Plot the data
    plt.plot(eta * beta_vals[beta] / K, UB, color='k')
    #plt.plot(eta, PhiB, color='r')

    # asymptotic solution eta < 2*K (phi) when beta = 1/2
    f_asym_in =  ( ( eta / K )**2 - 4.0 ) * (K / 4.0)
    #plt.plot(eta, f_asym_in, "g--")
    # asymptotic solution large eta
    f_asym_out = eta - 2*K
    #plt.plot(eta, f_asym_out, "g:")
    # uniformly valid asymptotic solution
    f_asym_uni = f_asym_in * 0.5 * ( np.tanh(-eta + 2.0*K) + 1) + f_asym_out * 0.5 * ( np.tanh(eta - 2.0*K) + 1)
    #plt.plot(eta, f_asym_uni, "b-")

    #plt.plot(eta, UBd, linestyle='dashed', color='k', label='UBd' )
    #plt.plot(eta, PsiB, color='g', label='PsiB' )
    #plt.plot(eta, ThetaB, color='r', label='ThetaB' )

    # Load numerical data
    zeta0s = [20, 40]
    lines = ["k--","k:","k-."]
    linecycler = cycle(lines)

    if plot_full_data :
        for i in range(0, len(zeta0s) ) :
            data2 = np.loadtxt("./DATA/K_" + str(K) + "_beta_" + str(beta_vals[beta]) + "_401x401_16_128" + "/Qout_" + str(zeta0s[i]) + ".dat")
            eta = data2[:,1]
            U = data2[:,8]
            eta_hzeta_0 = eta[0:401]
            U_hzeta_0 = U[0:401]
            plt.plot(eta_hzeta_0, U_hzeta_0, next(linecycler), label="U_num zeta0=" + str(zeta0s[i]) )

# axes etc
#plt.xlabel('eta')
axes = plt.gca()
#axes.set_xlim([0,40])
#axes.set_ylim([-K,10])
axes.set_xlim([0,1.5])
axes.set_ylim([0,1.2])

#axes.legend()

# Hide axis numbers
#axes.xaxis.set_ticklabels([])
#axes.yaxis.set_ticklabels([])

plt.savefig("profile_plot.eps", format='eps', dpi=1000)

plt.show()
