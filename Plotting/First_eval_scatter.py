import numpy as np
import matplotlib.pyplot as plt

K = 4
beta = 0

data = np.loadtxt("./DATA/K_"+ str(K) + "_zeta0_1_beta_"+ str(beta) + "_301x301_30_30/First_eval.dat")

alpha = data[:,0]
c_real = data[:,1]
c_imag = data[:,2]

plt.scatter(alpha, c_imag, c="black")
axes = plt.gca()
axes.set_xlim([0,0.7])
axes.set_ylim([0,0.1])
plt.grid(True)

plt.show()
