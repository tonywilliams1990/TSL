#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

zeta0 = 4
alpha = 1
beta = 0
K = 0

data = np.loadtxt("./DATA/K_"+ str(K) + "_alpha_" + str(alpha) + "_beta_"+ str(beta) + "_101x101_16_128/Qout_" + str(zeta0) + ".dat")

zeta_hat = data[:,0]
eta = data[:,1]
Phi = data[:,6]
Psi = data[:,7]
U = data[:,8]

zeta = zeta0 * zeta_hat
eta_hat = eta / zeta0

#V = (1/np.sqrt(2))*(eta*U - Phi)
#W = (1/np.sqrt(2))*(zeta_hat*U - Psi)

U[U > 1.0] = 1.0

min_x = np.min(zeta_hat)
max_x = np.max(zeta_hat)
min_y = np.min(eta)
max_y = np.max(eta)

npts = 500

xi = np.linspace(min_x, max_x, npts)
yi = np.linspace(min_y, max_y, npts)
Ui = mlab.griddata(zeta_hat, eta, U, xi, yi, interp = 'linear')
#Vi = mlab.griddata(zeta_hat, eta, V, xi, yi, interp = 'nn')
#Wi = mlab.griddata(zeta_hat, eta, W, xi, yi, interp = 'nn')

CS = plt.contourf(xi, yi, Ui, 10,
                  #[-1, -0.1, 0, 0.1],
                  #alpha=0.5,
                  cmap=plt.cm.YlGnBu_r)

CB = plt.colorbar(CS, shrink=0.9)

#plt.streamplot(xi, yi, Wi, Vi, arrowstyle='->', density=2.5, color='k')  

font_size = 22

plt.xlabel(r'$\hat{\zeta}$', fontsize=font_size)
plt.ylabel(r'$\eta$', rotation='horizontal', fontsize=font_size)
plt.title(r'$\beta =' + str(beta) + ', alpha=' + str(alpha) + ', \zeta_0=' + str(zeta0) + ', K=' + str(K) + '$', fontsize=font_size)


axes = plt.gca()
axes.set_xlim([0,2])
axes.set_ylim([0,40])

plt.show()
