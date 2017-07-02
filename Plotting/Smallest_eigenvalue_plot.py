#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle
import warnings
from scipy.interpolate import interp1d

plot_smallest_eigenvalue_data = False
plot_data_from_Rich = True
plot_critical_eigenvalue = True
hide_axis_labels = False

beta = 0.1

K_max = 3.0
K_min = -2.0
K_vals = np.linspace(K_min,K_max,51)

points = ["k^", "ks", "k*", "kv"]
pointcycler = cycle(points)

if plot_smallest_eigenvalue_data :

    # Plot eigenvalue (real part) for various beta
    beta_vals = [0.1, 0.2, 0.5, 0.6]
    for beta in beta_vals :
        point = next(pointcycler)
        for i in range( 0, len( K_vals ) ):
          # Load data
          with warnings.catch_warnings():
              warnings.simplefilter("ignore")
              data = np.loadtxt("./DATA/Eigenvalue_problem_beta_" + str(beta) + "/Eigenvalues_K_" + str(K_vals[i]) + ".dat")
          #print K_vals[i]
          eigenvals = data
          #print min_eigenval
          if eigenvals.size > 0:
              min_eigenval = np.min(eigenvals)
              plt.plot( K_vals[i], min_eigenval, point, clip_on=False )

    # Plot eigenvalue (real part) for beta = 0
    beta = 0
    K_vals_1 = np.arange(-2.0,0.0,0.1)
    K_vals_2 = np.arange(0.0,0.9,0.1)
    K_vals = np.concatenate((K_vals_1, K_vals_2), axis=0)
    point = "ko"
    for i in range( 0, len( K_vals ) ):
      # Load data
      data = np.loadtxt("./DATA/Eigenvalue_problem_beta_" + str(beta) + "/Eigenvalues_K_" + str(K_vals[i]) + ".dat")
      #print K_vals[i]
      eigenvals = data
      #print min_eigenval
      if eigenvals.size > 0:
          min_eigenval = np.min(eigenvals)
          plt.plot( K_vals[i], min_eigenval, point, clip_on=False )

    axes = plt.gca()
    axes.set_xlim([K_min,K_max])
    axes.set_ylim([-1,3])

    if hide_axis_labels :
        axes.xaxis.set_ticklabels([])
        axes.yaxis.set_ticklabels([])

    plt.savefig("Smallest_real_eigenvalues.eps", format='eps', dpi=1000)

if plot_data_from_Rich :
    plt.figure()
    beta_vals = [0.1, 0.2, 0.4, 0.5, 0.6]
    for beta in beta_vals :
        data = np.loadtxt("./DATA/Eigenvalue_data_from_Rich/evs_" + str(beta) + ".dat")
        K = -data[:,1]
        lambda_r = data[:,2]
        plt.plot( K, lambda_r, 'k')
        axes = plt.gca()
        axes.set_xlim([-2,3])
        axes.set_ylim([-1,3])
        if hide_axis_labels :
            axes.xaxis.set_ticklabels([])
            axes.yaxis.set_ticklabels([])
    plt.savefig("Rich_eigenvalue_data.eps", format='eps', dpi=1000)

if plot_critical_eigenvalue :
    plt.figure()
    K = [ -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6 ]
    #beta = [ -0.2917, -0.0187, 0.1944, 0.3473, 0.4419, 0.4862, 0.4990, 0.5009, 0.5011, 0.5012, 0.5013, 0.5014, 0.5015, 0.5016, 0.5017]
    beta = [ -0.2917, -0.0187, 0.1944, 0.3473, 0.4419, 0.4862, 0.4990, 0.4996, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    N_interp = 1000
    # Interpolate the data
    K_new = np.linspace( K[0], K[-1], N_interp, endpoint=True)
    beta_new = interp1d( K, beta, kind='cubic' )

    # Plot data
    #plt.plot( K, beta, 'ko')
    plt.plot( K_new, beta_new(K_new), 'k' )

    axes = plt.gca()
    axes.set_xlim([-1,6])
    axes.set_ylim([-0.3,0.6])

    if hide_axis_labels :
        axes.xaxis.set_ticklabels([])
        axes.yaxis.set_ticklabels([])

    plt.savefig("Critcal_eigenvalue.eps", format='eps', dpi=1000)

plt.show()
