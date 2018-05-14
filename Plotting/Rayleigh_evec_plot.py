#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


zeta0 = 1
beta = 0
K = 4
alpha = 0.53
number_of_levels = 11

data = np.loadtxt("./DATA/K_"+ str(K) + "_zeta0_" + str(zeta0) + "_beta_"+ str(beta) + "_201x201_30_30/alpha_"+ str(alpha) + "_evecs.dat")

zeta_hat = data[:,0]
eta = data[:,1]
evec1_real = data[:,2]
evec1_imag = data[:,3]
#evec1_real = data[:,8]
#evec1_imag = data[:,9]
evec1_abs = np.sqrt( evec1_real * evec1_real + evec1_imag * evec1_imag )

#TODO what if there is more than one eigenvalue ???

min_x = np.min(zeta_hat)
max_x = np.max(zeta_hat)
min_y = np.min(eta)
max_y = np.max(eta)

npts = 500

xi = np.linspace(min_x, max_x, npts)
yi = np.linspace(min_y, max_y, npts)
evec1_real_i = mlab.griddata(zeta_hat, eta, evec1_real, xi, yi, interp = 'linear')
evec1_imag_i = mlab.griddata(zeta_hat, eta, evec1_imag, xi, yi, interp = 'linear')
evec1_abs_i  = mlab.griddata(zeta_hat, eta, evec1_abs , xi, yi, interp = 'linear')

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

#plt.xlabel(r'$\hat{\zeta}$', fontsize=font_size)
#plt.ylabel(r'$\eta$', rotation='horizontal', fontsize=font_size)
#plt.title(r'$\beta =' + str(beta) + ', \zeta_0=' + str(zeta0) + ', K=' + str(K) + '$', fontsize=font_size)

axes = plt.gca()
axes.set_xlim([0,max_x])
#plt.xticks(np.arange(0, max_x + 0.5, 0.5))
axes.set_ylim([0,max_y])

plt.figure(2)

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

plt.figure(3)

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

plt.show()