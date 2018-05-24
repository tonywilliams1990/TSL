import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

#K_vals = [1, 2, 3, 4, 5]
K_vals = [8, 9, 10, 11, 12]
N = 301
beta = 0.5
markers = [".", "^", "o", "*", "v", "+", "D"]
markercycler = cycle(markers)

back_data = False

for K in K_vals :

    data = np.loadtxt("./DATA/K_"+ str(K) + "_zeta0_1_beta_"+ str(beta) + "_" + str(N) + "x" + str(N) + "_32_32/First_eval.dat")

    alpha = data[:,0]
    c_real = data[:,1]
    c_imag = data[:,2]

    marker = next(markercycler)
    plt.scatter(alpha, c_imag, c="black", marker=marker, clip_on=False)

    if back_data:
        data_back = np.loadtxt("./DATA/K_"+ str(K) + "_zeta0_1_beta_"+ str(beta) + "_" + str(N) + "x" + str(N) + "_32_32/First_eval_back.dat")
        alpha_back = data_back[:,0]
        c_real_back = data_back[:,1]
        c_imag_back = data_back[:,2]

        plt.scatter(alpha_back, c_imag_back, c="black", marker=marker, clip_on=False)


axes = plt.gca()
axes.set_xlim([0,1.0])
axes.set_ylim([0,0.08])
plt.grid(True)

plt.show()
