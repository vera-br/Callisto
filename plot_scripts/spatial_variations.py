import numpy as np
import matplotlib.pyplot as plt


def position_evolution_column(B_vectors1, times, color):
    B_components = B_vectors1.transpose()
    times = times / 3600  # convert to hours
    params = {
        'axes.labelsize': 14,
        'font.size': 18,
        'legend.fontsize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.figsize': [10, 10]
    }
    plt.rcParams.update(params)

    Bx = B_components[0]
    By = B_components[1]
    Bz = B_components[2]
    Bmag = np.linalg.norm(B_components, axis=0)

    fig, axs = plt.subplots(4, 1, sharex=True)
    linewidth=1.

    ax1 = axs[0]
    ax1.plot(times, Bx - np.mean(Bx), color=color, linewidth=linewidth)
    ax1.set_ylabel("X (nT)")

    ax1.set_xlim(times[0], times[-1])
    ax1.grid()
    ax1.minorticks_on()

    ax2 = axs[1]
    ax2.plot(times, By - np.mean(By), color=color, linewidth=linewidth)
    ax2.set_ylabel("Y (nT)")
    ax2.set_xlim(times[0], times[-1])
    ax2.grid()
    ax2.minorticks_on()

    ax3 = axs[2]
    ax3.plot(times, Bz - np.mean(Bz), color=color, linewidth=linewidth)
    ax3.set_ylabel("Z (nT)")
    ax3.set_xlim(times[0], times[-1])
    ax3.grid()
    ax3.minorticks_on()

    ax4 = axs[3]
    ax4.plot(times, Bmag - np.mean(Bmag), color=color, linewidth=linewidth)

    ax4.set_ylabel(r"$\vert$R$\vert$ (nT)")
    ax4.set_xlabel("Time (hours)")
    ax4.set_xlim(times[0], times[-1])
    ax4.grid()
    ax4.minorticks_on()

    fig.tight_layout()

    return plt.show()