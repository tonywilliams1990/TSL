#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

hide_axis_labels = False
save_fig = False
show_plots = True

K = 8
beta = 0
x_vals = [0, 5, 10, 100]
lines = ["k-", "k--","k:","k-."]
#lines = ["k-"]
linecycler = cycle(lines)
hzeta_max=16

for x in x_vals :
    data = np.loadtxt("./DATA/Isolated_Marching_K_" + str(K) + "_beta_0_1000x401x401_100_16_128/Wall_shear_zeta0_1_x_" + str(x) + ".dat")
    hzeta  = data[:,0]
    shear = data[:,1]

    line = next(linecycler)

    plt.figure(1)
    plt.plot(hzeta, shear, line)



plt.figure(1)
axes = plt.gca()
axes.set_xlim([0,hzeta_max])
if hide_axis_labels:
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
if save_fig:
    plt.savefig('./figs/Wall_shear_K_' + str(K) + '.eps', format='eps', dpi=1000)

if show_plots:
    plt.show()
