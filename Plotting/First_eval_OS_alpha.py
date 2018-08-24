import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

save_fig = False
hide_axis_labels = False
legend = True
legend_text_colour = "red"

R = 5000 * 5000
alpha_vals = [0.2, 0.4, 0.6, 0.8]
N = 601
beta = 0.5
zeta0 = 1
markers = [".", "^", "o", "*", "v", "+", "D"]
markercycler = cycle(markers)

for alpha in alpha_vals :

    data = np.loadtxt("./DATA/Eval_vs_K_alpha_"+ str(alpha) + ".dat")

    K = data[:,0]
    c_real = data[:,1]
    c_imag = data[:,2]

    marker = next(markercycler)
    plt.figure(1)
    plt.scatter(K, c_imag, c="black", marker=marker, clip_on=False, label=" a = " + str(alpha) )

    plt.figure(2)
    plt.scatter(K, c_real, c="black", marker=marker, clip_on=False, label=" a = " + str(alpha) )


plt.figure(1)
axes = plt.gca()
axes.set_xlim([7,12])
axes.set_ylim([-0.125,0.075])
axes.set_axisbelow(True)
plt.grid(True)

plt.figure(2)
axes = plt.gca()
axes.set_xlim([7,12])
axes.set_ylim([0.75,0.81])
axes.set_axisbelow(True)
plt.grid(True)

if legend:
    plt.figure(1)
    leg = plt.legend(framealpha=1)
    for text in leg.get_texts():
        text.set_color(legend_text_colour)

    plt.figure(2)
    leg = plt.legend(framealpha=1)
    for text in leg.get_texts():
        text.set_color(legend_text_colour)

if hide_axis_labels:
    plt.figure(1)
    axes = plt.gca()
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
    plt.figure(2)
    axes = plt.gca()
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
if save_fig:
    plt.figure(1)
    plt.savefig("./figs/First_eval_OS_alpha_imag.eps", format='eps', dpi=1000)
    plt.figure(2)
    plt.savefig("./figs/First_eval_OS_alpha_real.eps", format='eps', dpi=1000)

plt.show()
