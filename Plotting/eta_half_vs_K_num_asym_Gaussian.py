#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt



# Load data
data5 = np.loadtxt("./DATA/K_Step_beta_0_zeta0_5_Gaussian/A_file.dat")
K_5  = data5[:,1]
A_5  = data5[:,2]
eta_half_5   = data5[:,4]
hat_rs_5 = eta_half_5 / 5
K_max_5 = np.max(K_5)

data10 = np.loadtxt("./DATA/K_Step_beta_0_zeta0_10_Gaussian/A_file.dat")
K_10  = data10[:,1]
A_10  = data10[:,2]
eta_half_10   = data10[:,4]
hat_rs_10 = eta_half_10 / 10
K_max_10 = np.max(K_10)

data15 = np.loadtxt("./DATA/K_Step_beta_0_zeta0_15_Gaussian/A_file.dat")
K_15  = data15[:,1]
A_15  = data15[:,2]
eta_half_15   = data15[:,4]
hat_rs_15 = eta_half_15 / 15
K_max_15 = np.max(K_15)

data20 = np.loadtxt("./DATA/K_Step_beta_0_zeta0_20_Gaussian/A_file.dat")
K_20  = data20[:,1]
A_20  = data20[:,2]
eta_half_20   = data20[:,4]
K_max_20 = np.max(K_20)
hat_rs_20 = eta_half_20 / 20

# Plot the numerical data
#plt.plot(K_5, hat_rs_5, "r--")
#plt.plot(K_10, hat_rs_10, "k--")
#plt.plot(K_15, hat_rs_15, "g:")
#plt.plot(K_20, hat_rs_20, "b:")

#plt.xlabel('K')
#axes = plt.gca()
#axes.set_xlim([0,3])
#axes.set_ylim([0,5])

# Plot the asymptotic predictions
F_minus_inf = -1.238
K_crit = 2.6

K = np.linspace( K_crit, K_max_10)
#r_s = -( K ) / ( F_minus_inf * np.sqrt(np.pi) )
#plt.plot(K, r_s, color='r')

plt.figure()

# Plot the numerical data
plt.plot(K_5, A_5/25, "g--")
plt.plot(K_10, A_10/100, "k--")
plt.plot(K_15, A_15/225, "b--")
plt.plot(K_20, A_20/400, "r--")
# get new data for zeta0 = 20, 10 ??

# Plot the asymptotic predictions
A_zeta0_2 = - ( K * K ) / ( F_minus_inf * F_minus_inf * np.pi)
plt.plot(K, A_zeta0_2, color='r')

plt.xlabel('K')
plt.ylabel('A/zeta0^2')
axes = plt.gca()
axes.set_xlim([2.5,3])
#axes.set_ylim([-2,0])

#plt.savefig("plot.eps", format='eps', dpi=1000)

plt.show()
