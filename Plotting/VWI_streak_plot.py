#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


zeta0 = 1
beta = 0.5
N = 80
M = 80
K = 10.5
R = 25000000
Sigma = 0


#data = np.loadtxt("./DATA/VWI/K_" + str(K) + "_R_" + str(R) + "_Sigma_" + str(Sigma) + "_" + str(N+1) + "x" + str(M+1) + "_20_20.dat")
data = np.loadtxt("./DATA/VWI/K_" + str(K) + "_R_2.5e+07_Sigma_" + str(Sigma) + "_" + str(N+1) + "x" + str(M+1) + "_20_20.dat")

zeta_hat = data[:,0]
eta = data[:,1]
Phi_pert = data[:,10]
Psi_pert = data[:,11]
U_pert = data[:,12]
Theta_pert = data[:,13]
Phi = data[:,14]# Phi_real
Psi = data[:,15]# Psi_real
U = data[:,16]# U_real
Theta = data[:,17]# Theta_real

zeta = zeta0 * zeta_hat
eta_hat = eta / zeta0

V = (1 - beta)*eta*U - Phi
W = (1 - beta)*zeta_hat*U - Psi # effectively W / zeta0


#U[U > 1.0] = 1.0

min_x = np.min(zeta_hat)
#max_x = np.max(zeta_hat)
max_x = 10
min_y = np.min(eta)
#max_y = np.max(eta)
max_y = 10

npts = 500

xi = np.linspace(min_x, max_x, npts)
yi = np.linspace(min_y, max_y, npts)
Ui = mlab.griddata(zeta_hat, eta, U, xi, yi, interp = 'linear')
Vi = mlab.griddata(zeta_hat, eta, V, xi, yi, interp = 'linear')
Wi = mlab.griddata(zeta_hat, eta, W, xi, yi, interp = 'linear')
Vorticity_perturbation = mlab.griddata(zeta_hat, eta, Theta_pert * zeta0, xi, yi, interp = 'linear')

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

#CS2 = plt.contour(CS, levels=CS.levels[::1],
#                  colors='w',
#                  origin='lower',
#                  hold='on')

CB = plt.colorbar(CS, shrink=1)
CB.set_ticks(levels)

font_size = 22

#plt.xlabel(r'$\hat{\zeta}$', fontsize=font_size)
#plt.ylabel(r'$\eta$', rotation='horizontal', fontsize=font_size)
#plt.title(r'$\beta =' + str(beta) + ', \zeta_0=' + str(zeta0) + ', K=' + str(K) + '$', fontsize=font_size)

axes = plt.gca()
axes.set_xlim([0,max_x])
#plt.xticks(np.arange(0, max_x + 0.5, 0.5))
axes.set_ylim([0,max_y])


#plt.savefig('U_contour_K_' + str(K) + "_beta_" + str(beta) + "_zeta0_" + str(zeta0) + "_N_" + str(N) + ".eps", format='eps', dpi=1000)

plt.figure()

levels = np.linspace(-7, 0, 11)

CS = plt.contourf(xi, yi, Vorticity_perturbation, levels,
                  #[-1, -0.1, 0, 0.1],
                  #alpha=0.5,
                  cmap=cmap,
                  extend='both')


CB = plt.colorbar(CS, shrink=1)
CB.set_ticks(levels)

#speed = np.sqrt( Vi*Vi + Wi*Wi)
#lw = 5 * speed / speed.max()
lw = 1
# set linewidth = lw if you want linewidth varying with speed

plt.streamplot(xi, yi, Wi, Vi, arrowstyle='->', density=1, color='k', linewidth=lw)

#plt.xlabel(r'$\hat{\zeta}$', fontsize=font_size)
#plt.ylabel(r'$\eta$', rotation='horizontal', fontsize=font_size)
#plt.title(r'$\beta =' + str(beta) + ', \zeta_0=' + str(zeta0) + ', K=' + str(K) + '$', fontsize=font_size)

axes = plt.gca()
axes.set_xlim([0,max_x])
#plt.xticks(np.arange(0, max_x + 0.5, 0.5))
axes.set_ylim([0,max_y])


#plt.savefig('VW_stream_K_' + str(K) + "_beta_" + str(beta) + "_zeta0_" + str(zeta0) + "_N_" + str(N) + ".eps", format='eps', dpi=1000)

plt.show()
