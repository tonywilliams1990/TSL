#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle
import os, sys, errno

show_figures = False
beta = 0
#K_vals = [0,0.87,1.95]
K_vals = np.arange(0.0, 3.01, 0.01)

zeta0 = 20

#K_vals_test = np.arange(0.00, 3.01, 0.01)
#for i in range( 0, len( K_vals_test ) ):
#    print str(K_vals_test[i]).rstrip('0').rstrip('.')

#print K_vals_test

# Make a directory for storing the plots in
path = "./K_step_plot_beta_" + str(beta) + "_zeta0_" + str(zeta0)
try:
    os.mkdir( path, 0755 )
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

# Iterate over K values
for K in K_vals :

    # Load data
    data = np.loadtxt("./DATA/K_Step_beta_" + str(beta) + "_zeta0_" + str(zeta0) + "/Qout_K_" + str(K).rstrip('0').rstrip('.') + "_zeta0_" + str(zeta0) + ".dat")
    zeta_hat = data[:,0]
    eta = data[:,1]
    Phi_pert = data[:,2]
    Psi_pert = data[:,3]
    U_pert = data[:,4]
    Theta_pert = data[:,5]
    Phi = data[:,6]
    Psi = data[:,7] # Psi / zeta0
    U = data[:,8]
    Theta = data[:,9] # Theta / zeta0

    zeta = zeta0 * zeta_hat
    eta_hat = eta / zeta0

    V = (1 - beta)*eta*U - Phi
    W = (1 - beta)*zeta_hat*U - Psi # effectively W / zeta0

    U[U > 1.0] = 1.0

    min_x = np.min(zeta_hat)
    #max_x = np.max(zeta_hat)
    max_x = 2
    min_y = np.min(eta)
    #max_y = np.max(eta)
    max_y = max_x * zeta0

    npts = 500

    xi = np.linspace(min_x, max_x, npts)
    yi = np.linspace(min_y, max_y, npts)
    Ui = mlab.griddata(zeta_hat, eta, U, xi, yi, interp = 'linear')
    #Vi = mlab.griddata(zeta_hat, eta, V, xi, yi, interp = 'linear')
    #Wi = mlab.griddata(zeta_hat, eta, W, xi, yi, interp = 'linear')
    #Vorticity_perturbation = mlab.griddata(zeta_hat, eta, Theta_pert * zeta0, xi, yi, interp = 'linear')
    plt.figure()

    origin = 'lower'
    cmap = plt.cm.YlGnBu_r

    #levels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    levels = np.linspace(0,1,11)

    CS = plt.contourf(xi, yi, Ui, levels,
                      #[-1, -0.1, 0, 0.1],
                      #alpha=0.5,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    font_size = 22

    #plt.xlabel(r'$\hat{\zeta}$', fontsize=font_size)
    #plt.ylabel(r'$\eta$', rotation='horizontal', fontsize=font_size)
    #plt.title(r'$\beta =' + str(beta) + ', \zeta_0=' + str(zeta0) + ', K=' + str(K) + '$', fontsize=font_size)

    axes = plt.gca()
    axes.set_xlim([0,max_x])
    plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,40])

    plt.savefig(path + '/U_contour_K_' + str(K).rstrip('0').rstrip('.') + ".eps", format='eps', dpi=1000)
    print str(K).rstrip('0').rstrip('.')


# Hide axis numbers
#axes.xaxis.set_ticklabels([])
#axes.yaxis.set_ticklabels([])

if show_figures:
    plt.show()
