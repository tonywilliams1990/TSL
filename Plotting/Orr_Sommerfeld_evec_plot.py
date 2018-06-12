#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

show_fig_1 = False
show_fig_2 = False
show_fig_3 = False
show_vel = True

save_fig = True

zeta0 = 1
beta = 0.5
K = 9
R = 20000
alpha = 0.23
number_of_levels = 11

data = np.loadtxt("./DATA/K_"+ str(K) + "_zeta0_" + str(zeta0) + "_beta_"+ str(beta) + "_601x601_32_32/eigenvectors_R_"+ str(R*R) +"/alpha_"+ str(alpha) + "_evecs.dat")

zeta_hat = data[:,0]
eta = data[:,1]
evec1_real = data[:,2]# v_real
evec1_imag = data[:,3]# v_imag
evec2_real = data[:,4]# w_real
evec2_imag = data[:,5]# w_imag

evec1_abs = np.sqrt( evec1_real * evec1_real + evec1_imag * evec1_imag ) # abs(v)
evec2_abs = np.sqrt( evec2_real * evec2_real + evec2_imag * evec2_imag ) # abs(w)

v = evec1_abs / np.max(evec1_abs)
w = evec2_abs / np.max(evec2_abs)

vel = np.sqrt( v * v + w * w )
vel = vel / np.max(vel)

#TODO not sure if this is the correct in-plane velocity

#TODO what if there is more than one eigenvalue ???

min_x = np.min(zeta_hat)
#max_x = np.max(zeta_hat)
max_x = 8
min_y = np.min(eta)
#max_y = np.max(eta)
max_y = 8

npts = 500

xi = np.linspace(min_x, max_x, npts)
yi = np.linspace(min_y, max_y, npts)
evec1_real_i = mlab.griddata(zeta_hat, eta, evec1_real, xi, yi, interp = 'linear')
evec1_imag_i = mlab.griddata(zeta_hat, eta, evec1_imag, xi, yi, interp = 'linear')
evec1_abs_i  = mlab.griddata(zeta_hat, eta, evec1_abs , xi, yi, interp = 'linear')
vel_i        = mlab.griddata(zeta_hat, eta, vel       , xi, yi, interp = 'linear')

if show_fig_1:

    plt.figure(1)

    origin = 'lower'
    cmap = plt.cm.YlGnBu_r
    #levels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    levels = np.linspace(np.min(evec1_real),np.max(evec1_real),number_of_levels)

    CS = plt.contourf(xi, yi, evec1_real_i, levels,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    font_size = 22
    axes = plt.gca()
    axes.set_xlim([0,max_x])
    #plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,max_y])

if show_fig_2:

    plt.figure(2)
    origin = 'lower'
    cmap = plt.cm.YlGnBu_r
    levels = np.linspace(np.min(evec1_imag),np.max(evec1_imag),number_of_levels)

    CS = plt.contourf(xi, yi, evec1_imag_i, levels,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    axes = plt.gca()
    axes.set_xlim([0,max_x])
    #plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,max_y])

if show_fig_3:

    plt.figure(3)
    origin = 'lower'
    cmap = plt.cm.YlGnBu_r
    levels = np.linspace(np.min(evec1_abs),np.max(evec1_abs),number_of_levels)

    CS = plt.contourf(xi, yi, evec1_abs_i, levels,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    axes = plt.gca()
    axes.set_xlim([0,max_x])
    #plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,max_y])

if show_vel:

    plt.figure(4)
    origin = 'lower'
    cmap = plt.cm.YlOrRd
    levels = np.linspace(np.min(vel),np.max(vel),number_of_levels)

    CS = plt.contourf(xi, yi, vel_i, levels,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    axes = plt.gca()
    axes.set_xlim([0,max_x])
    #plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,max_y])

    if save_fig:
        plt.savefig("./figs/OS_2D_plane_vel_beta_" + str(beta) + "_R_" + str(R) + ".eps", format='eps', dpi=1000)

plt.show()
