import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from itertools import cycle

save_fig = True
hide_axis_labels = True
legend = True
legend_text_colour = "white"

#K_vals = [1, 2, 3, 4, 5]
K_vals = [8, 9, 10, 11, 12]
N = 601
beta = 0.5
zeta0 = 1
markers = [".", "^", "o", "*", "v", "+", "D"]
markercycler = cycle(markers)
tol = 0.003

back_data = False

alpha_max_growth = 0
max_growth_rate = 0

print "K\talpha_max\tmax_growth_rate"

for K in K_vals :

    data = np.loadtxt("./DATA/Inviscid_stability_601_DATA/K_"+ str(K) + "_zeta0_1_beta_"+ str(beta) + "_" + str(N) + "x" + str(N) + "_32_32/First_eval.dat")

    #alpha = data[:,0]
    #c_real = data[:,1]
    #c_imag = data[:,2]
    # Filter out c_i < tol
    c_imag_filter = np.abs(data[:,2])>tol
    alpha = data[:,0][c_imag_filter]
    c_real = data[:,1][c_imag_filter]
    c_imag = data[:,2][c_imag_filter]

    max_growth_rate = np.max( alpha * c_imag )
    index = np.argmax( alpha * c_imag )
    alpha_max_growth = alpha[index]

    marker = next(markercycler)
    plt.scatter(alpha, c_imag, c="black", marker=marker, clip_on=False, label=" K = " + str(K) )

    if back_data:
        data_back = np.loadtxt("./DATA/Inviscid_stability_601_DATA/K_"+ str(K) + "_zeta0_1_beta_"+ str(beta) + "_" + str(N) + "x" + str(N) + "_32_32/First_eval_back.dat")
        c_imag_back_filter = np.abs(data_back[:,2])>tol
        alpha_back = data_back[:,0][c_imag_back_filter]
        c_real_back = data_back[:,1][c_imag_back_filter]
        c_imag_back = data_back[:,2][c_imag_back_filter]

        plt.scatter(alpha_back, c_imag_back, c="black", marker=marker, clip_on=False)

        if ( np.max( alpha_back * c_imag_back ) > max_growth_rate):
            max_growth_rate = np.max( alpha_back * c_imag_back )
            index = np.argmax( alpha_back * c_imag_back )
            alpha_max_growth = alpha_back[index]

    print K, "\t", alpha_max_growth, "\t\t", max_growth_rate

axes = plt.gca()
#axes.set_xlim([0,0.65])
#axes.set_ylim([0,0.12])
axes.set_xlim([0,1.0])
axes.set_ylim([0,0.07])
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
    plt.savefig("./figs/First_eval_scatter_beta_" + str(beta) + "_zeta0_" + str(zeta0) + ".eps", format='eps', dpi=1000)

plt.show()
