#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

hide_axis_labels = True
save_fig = True
show_plots = False
plot_lines = False
beta=0
zeta0 = 1
K_vals = [1, 2]
lines = ["k--","k:","k-.", "k-"] #solid black (preffered for thesis)
linecycler = cycle(lines)


x_max = 100

for K in K_vals :
    data = np.loadtxt("./DATA/Marching_K_" + str(K) + "_zeta0_" + str(zeta0) + "_2000x401x401_100_16_128/A_file.dat")

    A = data[:,1]
    x = data[:,2]
    shear_at_origin = data[:,3]
    streak_radius = data[:,4]
    integral_U_squared = data[:,5]

    line = next(linecycler)

    # A plot
    plt.figure(1)
    plt.plot(x, A, line)

    # Shear plot
    plt.figure(2)
    plt.plot(x, shear_at_origin, line)

    # integral U^2 plot
    plt.figure(3)
    plt.plot(x, integral_U_squared, line)

    data2 = np.loadtxt("./DATA/Marching_K_" + str(K) + "_zeta0_" + str(zeta0) + "_1000x201x201_100_16_64/A_file.dat")

    A_2 = data2[:,1]
    x_2 = data2[:,2]
    shear_at_origin_2 = data2[:,3]
    streak_radius_2 = data2[:,4]
    integral_U_squared_2 = data2[:,5]

    line = next(linecycler)

    # A plot
    plt.figure(1)
    plt.plot(x_2, A_2, line)

    # Shear plot
    plt.figure(2)
    plt.plot(x_2, shear_at_origin_2, line)

    # integral U^2 plot
    plt.figure(3)
    plt.plot(x_2, integral_U_squared_2, line)

# Set axes limits
plt.figure(1)
axes = plt.gca()
axes.set_xlim([0,x_max])
if hide_axis_labels:
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
if save_fig:
    plt.savefig('./figs/A_vs_x_zeta0_' + str(zeta0) + '.eps', format='eps', dpi=1000)
plt.figure(2)
axes = plt.gca()
axes.set_xlim([0,x_max])
if hide_axis_labels:
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
if save_fig:
    plt.savefig('./figs/shear_at_origin_vs_x_zeta0_' + str(zeta0) + '.eps', format='eps', dpi=1000)
plt.figure(3)
axes = plt.gca()
axes.set_xlim([0,x_max])
if hide_axis_labels:
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
if save_fig:
    plt.savefig('./figs/integral_U_squared_vs_x_zeta0_' + str(zeta0) + '.eps', format='eps', dpi=1000)


if show_plots:
    plt.show()
