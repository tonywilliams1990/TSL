#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

Hide_axis = True

zeta0 = 20
beta = 0
K_val = [0.2, 1.2]
N = 1

for j in range(0, len(K_val) ) :

    data = np.loadtxt("./DATA/K_Step_beta_"+ str(beta) + "_zeta0_" + str(zeta0) + "_N_" + str(N) + "/Qout_K_" + str(K_val[j]) + "_zeta0_" + str(zeta0) + ".dat")

    zeta_hat = data[:,0]
    eta = data[:,1]
    Phi_pert = data[:,2]
    Psi_pert = data[:,3]
    U_pert = data[:,4]
    Theta_pert = data[:,5]
    Phi = data[:,6]
    Psi = data[:,7]
    U = data[:,8]
    Theta = data[:,9]

    index = 0
    for i in range( 0, len(zeta_hat) ):
        if zeta_hat[i] < 0.7 :
            index = index + 1

    #print index

    zeta = zeta0 * zeta_hat
    vals = 400
    Psi_scaled = Psi[index: index + vals] / 0.7
    zeta_scaled = zeta[index: index + vals]
    eta_scaled = eta[index: index + vals]

    plt.plot(eta_scaled, Psi_scaled, "k-", linewidth=1)

axes = plt.gca()
axes.set_xlim([0,zeta0])
axes.set_ylim([0,1.2])

if Hide_axis :
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])

plt.savefig("Cross_flow.eps", format='eps', dpi=1000)

plt.show()
