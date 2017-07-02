#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

hide_axis_labels = False

zeta0 = 20
beta = 0
N = 400
K_vals = [0.9, 1.4, 1.9]
lines = [ 'solid', 'dashdot', 'dotted', 'dashed' ]
linecycler = cycle(lines)

for i in range( 0, len( K_vals ) ) :
    K = K_vals[i]

    data = np.loadtxt("./DATA/K_Step_beta_"+ str(beta) + "_zeta0_" + str(zeta0) + "/Qout_K_" + str(K) + "_zeta0_" + str(zeta0) + ".dat")

    zeta_hat = data[:,0]
    eta = data[:,1]
    U = data[:,8]

    zeta = zeta0 * zeta_hat
    eta_hat = eta / zeta0

    # Find value of eta on zeta=0 at which U=1/2
    eta_centreline = eta[0:N+1]
    U_centreline = U[0:N+1]
    # find the nodes either side of U=0.5
    lower = 0
    upper = 1
    for j in range( 0, len( eta_centreline ) ) :
        if U_centreline[j] < 0.5 and U_centreline[j+1] > 0.5 :
            lower = j
            upper = j+1

    # linearly interpolate
    r_s =  ( 0.5 - U_centreline[lower] ) * ( eta_centreline[upper] - eta_centreline[lower] ) / ( U_centreline[upper] - U_centreline[lower]  ) + eta_centreline[lower]
    print r_s
    zeta_r_s = zeta / r_s
    eta_r_s = eta / r_s

    min_x = np.min(zeta_r_s)
    #max_x = np.max(zeta_hat)
    max_x = 1.5
    min_y = np.min(eta_r_s)
    #max_y = np.max(eta)
    max_y = 1.5

    npts = 500

    xi = np.linspace(min_x, max_x, npts)
    yi = np.linspace(min_y, max_y, npts)
    Ui = mlab.griddata(zeta_r_s, eta_r_s, U, xi, yi, interp = 'linear')

    levels = [ 0.5 ]
    #levels = np.linspace(0,1,11)

    CS = plt.contour(xi, yi, Ui, levels,
                            colors='k',
                            linestyles=next(linecycler),
                            origin='lower',
                            extend='both')

# axes
axes = plt.gca()
axes.set_xlim([0,1.2])
axes.set_ylim([0,1.2])

if hide_axis_labels :
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])

# plot a circle idicating predicted circular shear layer for strong injection
circle = plt.Circle((0,0), 1, color='k', linestyle=next(linecycler), fill=False)
axes.add_artist(circle)

plt.savefig( "Moderate_injection_circles.eps", format='eps', dpi=1000 )

plt.show()
