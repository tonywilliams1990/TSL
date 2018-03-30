#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

hide_axis_labels = False
save_fig = False
show_plots = True
plot_lines = False
beta=0
#K = 2
K_vals = [0.5, 1, 2]
zeta0_vals = [1, 2]
#lines = ["k--","k:","k-.", "k-"] #solid black (preffered for thesis)
lines = ["k--", "k-"]
lines2 = ["k--","k:","k-."] # dashed, dotted, dash-dotted
linecycler = cycle(lines)
#lines2 = ["b-.","r-.","g-.","y-."]
linecycler2 = cycle(lines2)

x_max = 100

for K in K_vals :
    for zeta0 in zeta0_vals :

        #data = np.loadtxt("./DATA/Marching_K_"+ str(K) + "_beta_" + str(beta) + "_1000x201x201_100_30_30/A_file.dat")
        #data = np.loadtxt("./DATA/Marching_DATA_1000x401x401_100_16_128/Marching_K_"+ str(K) + "_zeta0_" + str(zeta0) + "/A_file.dat")
        #data = np.loadtxt("./DATA/Isolated_Marching_K_" + str(K) + "_beta_0_1000x201x201_100_16_64/A_file.dat")
        #data = np.loadtxt("./DATA/Isolated_Marching_K_" + str(K) + "_beta_0_1000x401x401_100_16_128/A_file.dat")
        data = np.loadtxt("./DATA/Marching_K_" + str(K) + "_zeta0_" + str(zeta0) + "_1000x201x201_100_16_64/A_file.dat")

        #zeta0 = data[:,0]
        A = data[:,1]
        x = data[:,2]
        shear_at_origin = data[:,3]
        streak_radius = data[:,4]
        integral_U_squared = data[:,5]

        if plot_lines:
            data2 = np.loadtxt("./DATA/K_" + str(K) + "_beta_" + str(beta) + "_201x201_16_64/A_file.dat")
            zeta0_sim = data2[0]
            A_sim = data2[1]
            U_eta_sim = data2[2]
            int_U_2_sim = data2[4]

        line = next(linecycler)
        line2 = next(linecycler2)

        # A plot
        plt.figure(1)
        plt.plot(x, A, line)
        if plot_lines:
            plt.plot([0, x_max], [A_sim, A_sim], line2)

        # Shear plot
        plt.figure(2)
        plt.plot(x, shear_at_origin, line)
        if plot_lines:
            plt.plot([0, x_max], [U_eta_sim, U_eta_sim], line2)

        # integral U^2 plot
        plt.figure(3)
        plt.plot(x, integral_U_squared, line)
        if plot_lines:
            plt.plot([0, x_max], [int_U_2_sim, int_U_2_sim], line2)

# Set axes limits
plt.figure(1)
axes = plt.gca()
axes.set_xlim([0,x_max])
if hide_axis_labels:
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
if save_fig:
    plt.savefig('./figs/A_vs_x.eps', format='eps', dpi=1000)
plt.figure(2)
axes = plt.gca()
axes.set_xlim([0,x_max])
if hide_axis_labels:
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
if save_fig:
    plt.savefig('./figs/shear_at_origin_vs_x.eps', format='eps', dpi=1000)
plt.figure(3)
axes = plt.gca()
axes.set_xlim([0,x_max])
if hide_axis_labels:
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
if save_fig:
    plt.savefig('./figs/integral_U_squared_vs_x.eps', format='eps', dpi=1000)


if show_plots:
    plt.show()
