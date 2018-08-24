#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import re

show_fig_1 = False
show_fig_2 = False
show_p_abs = False
show_u_abs = True
show_v_abs = False
show_w_abs = False
show_vel = True

save_fig = True

N = 600
N_p_1 = N + 1
zeta0 = 1
beta = 0.5
K = 9
alpha = 0.23
number_of_levels = 11


data = np.loadtxt("./DATA/Inviscid_stability_601_DATA/K_"+ str(K) + "_zeta0_" + str(zeta0) + "_beta_"+ str(beta) + "_" + str(N_p_1) + "x" + str(N_p_1) + "_32_32/eigenvectors/alpha_"+ str(alpha) + "_evecs.dat")
streak_data = np.loadtxt("./DATA/Inviscid_stability_601_DATA/K_"+ str(K) + "_zeta0_" + str(zeta0) + "_beta_"+ str(beta) + "_" + str(N_p_1) + "x" + str(N_p_1) + "_32_32/Qout_"+ str(zeta0) +".dat")
# Get the eigenvalue
fname = "./DATA/Inviscid_stability_601_DATA/K_"+ str(K) + "_zeta0_" + str(zeta0) + "_beta_"+ str(beta) + "_" + str(N_p_1) + "x" + str(N_p_1) + "_32_32/eigenvalues/alpha_"+ str(alpha) + "_evals.dat"
eigenvalue_data = np.loadtxt(fname, dtype = str)
eval_string = eigenvalue_data.tostring()
eval_string = eval_string[1:]
eval_string = eval_string[:-1]
c_r = float(eval_string.split(',')[0])
c_i = float(eval_string.split(',')[1])
c = complex(c_r, c_i)

zeta_hat = data[:,0]
eta = data[:,1]
evec1_real = data[:,2]
evec1_imag = data[:,3]
U = streak_data[:,8]
evec1_abs = np.sqrt( evec1_real * evec1_real + evec1_imag * evec1_imag )
evec1_max = np.max(evec1_abs)
evec1_real = evec1_real / evec1_max
evec1_imag = evec1_imag / evec1_max

p_real_eta = np.zeros(len(evec1_real))
p_imag_eta = np.zeros(len(evec1_real))
p_eta_abs_2 = np.zeros(len(evec1_real))
U_eta = np.zeros(len(U))

p_real_zeta = np.zeros(len(evec1_real))
p_imag_zeta = np.zeros(len(evec1_real))
p_zeta_abs_2 = np.zeros(len(evec1_real))
U_zeta = np.zeros(len(U))

min_d_zeta = 1.0

for i in range(N):
    for j in range(N):
        if( i != 0 and j !=0 and i != N and j != N ):
            d_eta = eta[i * N + j + 1] - eta[i * N + j - 1]
            p_real_eta[i * N + j] = (evec1_real[i * N + j + 1] - evec1_real[i * N + j - 1]) / (2 * d_eta)
            p_imag_eta[i * N + j] = (evec1_imag[i * N + j + 1] - evec1_imag[i * N + j - 1]) / (2 * d_eta)
            p_eta_abs_2[i * N + j] = p_real_eta[i * N + j] * p_real_eta[i * N + j] + p_imag_eta[i * N + j] * p_imag_eta[i * N + j]
            U_eta[i * N + j] = (U[i * N + j + 1] - U[i * N + j - 1]) / (2 * d_eta)

for i in range(N):
    for j in range(N):
        if( i != 0 and j !=0 and i != N and j != N ):
            d_zeta = zeta_hat[(i + 1) * N + j] - zeta_hat[(i - 1) * N + j]
            p_real_zeta[i * N + j] = (evec1_real[(i + 1) * N + j] - evec1_real[(i - 1) * N + j]) / (2 * d_zeta)
            p_imag_zeta[i * N + j] = (evec1_imag[(i + 1) * N + j] - evec1_imag[(i - 1) * N + j]) / (2 * d_zeta)
            p_zeta_abs_2[i * N + j] = p_real_zeta[i * N + j] * p_real_zeta[i * N + j] + p_imag_zeta[i * N + j] * p_imag_zeta[i * N + j]
            U_zeta[i * N + j] = ( U[(i + 1) * N + j] - U[(i - 1) * N + j] ) / (2 * d_zeta)


p = evec1_real + 1j * evec1_imag
p_eta = p_real_eta + 1j * p_imag_eta
p_zeta = p_real_zeta + 1j * p_imag_zeta


U_zeta = U_zeta / np.max(U_zeta)
U_eta = U_eta / np.max(U_eta)

div = 1j * alpha * (U - c)
v = - p_eta / div
w = - p_zeta / div

u = - 1j * alpha * p - U_eta * v - U_zeta * w
u = u / ( 1j * alpha * (U - c) )

#rhs = - 1j * alpha * p - U_eta * v - U_zeta * w
#rhs_conj = 1j * alpha * np.conjugate(p) - U_eta * np.conjugate(v) - U_zeta * np.conjugate(w)
#u_abs = rhs * rhs_conj
#u_abs = u_abs / ( alpha * alpha * ( U - c ) * ( U - np.conjugate(c) ) )

p_abs = evec1_abs / np.max(evec1_abs)

U_c_conj = U * U
U_c_conj = U_c_conj - 2.0 * c_r * U
U_c_conj = U_c_conj + c_r * c_r + c_i * c_i
U_c_conj = alpha * alpha * U_c_conj


#u_abs = evec1_abs / U_c_conj
#u_abs = np.real(u_abs) / U_c_conj
#u_abs = np.sqrt(u_abs)
u_abs = np.absolute(u)
u_abs = u_abs / np.max(u_abs)

#v_abs = evec1_abs / U_c_conj
#v_abs = p_eta_abs_2 / U_c_conj
v_abs = np.absolute(v)
#v_abs = np.sqrt(np.real(v_abs))
v_abs = v_abs / np.max(v_abs)

#w_abs = p_zeta_abs_2 / U_c_conj
w_abs = np.absolute(w)
w_abs = w_abs / np.max(w_abs)

vel = np.sqrt( v_abs * v_abs + w_abs * w_abs )
#vel = np.absolute( np.sqrt( v * v + w * w ) )
vel = vel / np.max(vel)

min_x = np.min(zeta_hat)
max_x = np.max(zeta_hat)
min_y = np.min(eta)
max_y = np.max(eta)

npts = 1000

xi = np.linspace(min_x, max_x, npts)
yi = np.linspace(min_y, max_y, npts)


if show_fig_1:
    evec1_real_i = mlab.griddata(zeta_hat, eta, evec1_real, xi, yi, interp = 'linear')
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
    evec1_imag_i = mlab.griddata(zeta_hat, eta, evec1_imag, xi, yi, interp = 'linear')
    plt.figure(2)
    origin = 'lower'
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

if show_p_abs:
    p_i  = mlab.griddata(zeta_hat, eta, p_abs , xi, yi, interp = 'linear')
    plt.figure(3)
    origin = 'lower'
    cmap = plt.cm.YlGnBu

    levels = np.linspace(np.min(p_abs),np.max(p_abs),number_of_levels)

    CS = plt.contourf(xi, yi, p_i, levels,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    axes = plt.gca()
    axes.set_xlim([0,16])
    #plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,16])

if show_u_abs:
    u_i = mlab.griddata(zeta_hat, eta, u_abs, xi, yi, interp = 'linear')
    plt.figure(4)
    origin = 'lower'
    cmap = plt.cm.YlOrRd


    levels = np.linspace(np.min(u_abs),np.max(u_abs),11)

    CS = plt.contourf(xi, yi, u_i, levels,
                      #[-1, -0.1, 0, 0.1],
                      #alpha=0.5,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    axes = plt.gca()
    axes.set_xlim([0,6])
    #plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,6])

    if save_fig:
        plt.savefig("./figs/Rayleigh_2D_u_abs_beta_" + str(beta) + "_zeta0_" + str(zeta0) + ".eps", format='eps', dpi=1000)

if show_v_abs:
    v_i = mlab.griddata(zeta_hat, eta, v_abs, xi, yi, interp = 'linear')
    plt.figure(5)
    origin = 'lower'
    cmap = plt.cm.YlOrRd


    levels = np.linspace(0,1,11)

    CS = plt.contourf(xi, yi, v_i, levels,
                      #[-1, -0.1, 0, 0.1],
                      #alpha=0.5,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    axes = plt.gca()
    axes.set_xlim([0,15])
    #plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,15])

if show_w_abs:
    w_i = mlab.griddata(zeta_hat, eta, w_abs, xi, yi, interp = 'linear')
    plt.figure(6)
    origin = 'lower'
    cmap = plt.cm.YlOrRd


    levels = np.linspace(0,1,11)

    CS = plt.contourf(xi, yi, w_i, levels,
                      #[-1, -0.1, 0, 0.1],
                      #alpha=0.5,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    axes = plt.gca()
    axes.set_xlim([0,15])
    #plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,15])

if show_vel:
    vel_i = mlab.griddata(zeta_hat, eta, vel, xi, yi, interp = 'linear')
    plt.figure(7)
    origin = 'lower'
    cmap = plt.cm.YlOrRd


    levels = np.linspace(0,1,11)

    CS = plt.contourf(xi, yi, vel_i, levels,
                      #[-1, -0.1, 0, 0.1],
                      #alpha=0.5,
                      cmap=cmap,
                      origin=origin,
                      extend='both')

    CB = plt.colorbar(CS, shrink=1)
    CB.set_ticks(levels)

    axes = plt.gca()
    axes.set_xlim([0,6])
    #plt.xticks(np.arange(0, max_x + 0.5, 0.5))
    axes.set_ylim([0,6])

    if save_fig:
        plt.savefig("./figs/Rayleigh_2D_plane_vel_beta_" + str(beta) + "_zeta0_" + str(zeta0) + ".eps", format='eps', dpi=1000)

plt.show()
