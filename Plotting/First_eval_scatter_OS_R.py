import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

save_fig = True
hide_axis_labels = True
legend = True
legend_text_colour = "white"

K = 9
R_vals = [1000, 2500, 5000, 10000, 20000]
N = 601
beta = 0.5
zeta0 = 1
markers = [".", "^", "o", "*", "v", "+", "D"]
markercycler = cycle(markers)
alpha_min = 0.04

back_data = True

alpha_max_growth = 0
max_growth_rate = 0

print "K\talpha_max\tmax_growth_rate"

for R in R_vals :

    data = np.loadtxt("./DATA/Viscous_stability_601_DATA/K_"+ str(K) + "_zeta0_1_beta_"+ str(beta) + "_" + str(N) + "x" + str(N) + "_32_32/First_eval_R_" + str(R * R) + ".dat")

    alpha_filter = np.abs(data[:,0])>alpha_min
    alpha = data[:,0][alpha_filter]
    c_real = data[:,1][alpha_filter]
    c_imag = data[:,2][alpha_filter]

    max_growth_rate = np.max( alpha * c_imag )
    index = np.argmax( alpha * c_imag )
    alpha_max_growth = alpha[index]

    marker = next(markercycler)
    plt.scatter(alpha, c_imag, c="black", marker=marker, clip_on=False, label="Rx^1/2 = " + str(R) )

    if back_data:
        data_back = np.loadtxt("./DATA/Viscous_stability_601_DATA/K_"+ str(K) + "_zeta0_1_beta_"+ str(beta) + "_" + str(N) + "x" + str(N) + "_32_32/First_eval_back_R_" + str(R * R) + ".dat")
        alpha_back_filter = np.abs(data_back[:,0])>alpha_min
        alpha_back = data_back[:,0][alpha_back_filter]
        c_real_back = data_back[:,1][alpha_back_filter]
        c_imag_back = data_back[:,2][alpha_back_filter]

        alpha_back = alpha_back[0:8]
        c_imag_back = c_imag_back[0:8]

        plt.scatter(alpha_back, c_imag_back, c="black", marker=marker, clip_on=False)

        if ( np.max( alpha_back * c_imag_back ) > max_growth_rate):
            max_growth_rate = np.max( alpha_back * c_imag_back )
            index = np.argmax( alpha_back * c_imag_back )
            alpha_max_growth = alpha_back[index]

    print R, "\t", alpha_max_growth, "\t\t", max_growth_rate

axes = plt.gca()
axes.set_xlim([0,0.5])
axes.set_ylim([-0.02,0.03])
axes.set_axisbelow(True)
plt.grid(True)

if legend:
    leg = plt.legend(framealpha=1)
    for text in leg.get_texts():
        text.set_color(legend_text_colour)

if hide_axis_labels:
    axes.xaxis.set_ticklabels([])
    axes.yaxis.set_ticklabels([])
if save_fig:
    plt.savefig("./figs/First_eval_scatter_beta_" + str(beta) + "_K_" + str(K) + ".eps", format='eps', dpi=1000)

plt.show()
