#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

show_fig_1 = False
show_fig_2 = False
show_fig_3 = False
show_vel = True

save_fig = False

zeta0 = 1
beta = 0.5
number_of_levels = 11
N = 80
M = 80
K = 10.5
R = 25000000
Sigma = 0

#data = np.loadtxt("./DATA/VWI/K_" + str(K) + "_R_" + str(R) + "_Sigma_" + str(Sigma) + "_" + str(N+1) + "x" + str(M+1) + "_20_20.dat")
data = np.loadtxt("./DATA/VWI/K_" + str(K) + "_R_2.5e+07_Sigma_" + str(Sigma) + "_" + str(N+1) + "x" + str(M+1) + "_20_20.dat")

zeta_hat = data[:,0]
eta = data[:,1]
evec1_real = data[:,2]# v_real
evec1_imag = data[:,3]# v_imag
evec2_real = data[:,4]# w_real
evec2_imag = data[:,5]# w_imag
evec3_real = data[:,6]# q_real
evec3_imag = data[:,7]# q_imag
evec4_real = data[:,8]# s_real
evec4_imag = data[:,9]# s_imag

evec1_abs = np.sqrt( evec1_real * evec1_real + evec1_imag * evec1_imag ) # abs(v)
evec2_abs = np.sqrt( evec2_real * evec2_real + evec2_imag * evec2_imag ) # abs(w)
evec3_abs = np.sqrt( evec3_real * evec3_real + evec3_imag * evec3_imag ) # abs(q)
evec4_abs = np.sqrt( evec4_real * evec4_real + evec4_imag * evec4_imag ) # abs(s)

v = evec1_abs / np.max(evec1_abs)
w = evec2_abs / np.max(evec2_abs)

vel = np.sqrt( v * v + w * w )
vel = vel / np.max(vel)


min_x = np.min(zeta_hat)
#max_x = np.max(zeta_hat)
max_x = 10
min_y = np.min(eta)
#max_y = np.max(eta)
max_y = 10

npts = 500

xi = np.linspace(min_x, max_x, npts)
yi = np.linspace(min_y, max_y, npts)
evec1_real_i = mlab.griddata(zeta_hat, eta, evec1_real, xi, yi, interp = 'linear')
evec1_imag_i = mlab.griddata(zeta_hat, eta, evec1_imag, xi, yi, interp = 'linear')
evec1_abs_i  = mlab.griddata(zeta_hat, eta, evec1_abs , xi, yi, interp = 'linear')
evec2_real_i = mlab.griddata(zeta_hat, eta, evec2_real, xi, yi, interp = 'linear')
evec2_imag_i = mlab.griddata(zeta_hat, eta, evec2_imag, xi, yi, interp = 'linear')
evec2_abs_i  = mlab.griddata(zeta_hat, eta, evec2_abs , xi, yi, interp = 'linear')
evec3_real_i = mlab.griddata(zeta_hat, eta, evec3_real, xi, yi, interp = 'linear')
evec3_imag_i = mlab.griddata(zeta_hat, eta, evec3_imag, xi, yi, interp = 'linear')
evec3_abs_i  = mlab.griddata(zeta_hat, eta, evec3_abs , xi, yi, interp = 'linear')
evec4_real_i = mlab.griddata(zeta_hat, eta, evec4_real, xi, yi, interp = 'linear')
evec4_imag_i = mlab.griddata(zeta_hat, eta, evec4_imag, xi, yi, interp = 'linear')
evec4_abs_i  = mlab.griddata(zeta_hat, eta, evec4_abs , xi, yi, interp = 'linear')
vel_i        = mlab.griddata(zeta_hat, eta, vel       , xi, yi, interp = 'linear')

if show_fig_1:

    plt.figure(1)

    origin = 'lower'
    cmap = plt.cm.YlGnBu_r
    #levels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    levels = np.linspace(np.min(evec3_real),np.max(evec3_real),number_of_levels)

    CS = plt.contourf(xi, yi, evec3_real_i, levels,
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
    levels = np.linspace(np.min(evec3_imag),np.max(evec3_imag),number_of_levels)

    CS = plt.contourf(xi, yi, evec3_imag_i, levels,
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
    levels = np.linspace(np.min(evec3_abs),np.max(evec3_abs),number_of_levels)

    CS = plt.contourf(xi, yi, evec3_abs_i, levels,
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
