import matplotlib.pyplot as plt
import numpy as np


def plot_B_evolution_cylindrical(B_vectors, times):
    B_components = B_vectors.transpose()
    times = times / 3600  # convert to hours
    params = {
        "axes.labelsize": 14,
        "font.size": 18,
        "legend.fontsize": 22,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "figure.figsize": [14, 6],
    }
    plt.rcParams.update(params)

    Brho = B_components[0]
    Bphi = B_components[1]
    Bz = B_components[2]
    Bmag = np.sqrt(Brho**2 + Bphi**2 + Bz**2)

    fig, axs = plt.subplots(4, 1, sharex=True)
    ax1 = axs[0]
    ax1.plot(times, Brho)
    ax1.set_ylabel(r"B$\rho$ (nT)")

    ax2 = axs[1]
    ax2.plot(times, Bphi)
    ax2.set_ylabel(r"B$\phi$ (nT)")

    ax3 = axs[2]
    ax3.plot(times, Bz)
    ax3.set_ylabel(r"Bz (nT)")
    ax3.set_xlabel("time (hr)")

    ax4 = axs[3]
    ax4.plot(times, Bmag)
    ax4.set_ylabel(r"$\vert$B$\vert$ (nT)")
    ax4.set_xlabel("time (hr)")

    fig.tight_layout()

    return plt.show()


def plot_B_evolution_poster_ind(
    B_vectors1,
    times,
    B_vectors2: np.ndarray = None,
    B_vectors3: np.ndarray = None,
    label1: str = None,
    label2: str = None,
    label3: str = None,
    filename: str = None,
    color1: str = None,
    color2: str = None,
    color3: str = None,
    linestyle1: str = None,
    linestyle2: str = None,
    linestyle3: str = None,
    linewidth1: float = None,
    linewidth2: float = None,
    bbox_coords: tuple = None,
    pad: float = None,
):
    B_components = B_vectors1.transpose()
    times = times / 3600  # convert to hours
    params = {
        "axes.labelsize": 20,
        "font.size": 18,
        "legend.fontsize": 22,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "figure.figsize": [12, 6],
    }
    plt.rcParams.update(params)

    Bx = B_components[0]
    By = B_components[1]
    Bz = B_components[2]
    Bmag = np.linalg.norm(B_components, axis=0)

    fig, axs = plt.subplots(ncols=2, nrows=1, sharex=True)

    ax1 = axs[0]
    ax1.plot(
        times,
        Bx,
        label=label1,
        color=color1,
        linewidth=linewidth1,
        linestyle=linestyle1,
    )
    ax1.set_ylabel("Bx (nT)")
    ax1.set_xlim(times[0], times[-1])
    ax1.grid()
    ax1.minorticks_on()

    ax2 = axs[1]
    ax2.plot(times, By, color=color1, linewidth=linewidth1, linestyle=linestyle1)
    ax2.set_ylabel("By (nT)")
    ax2.set_xlim(times[0], times[-1])
    ax2.grid()
    ax2.minorticks_on()

    alpha = 1.
    if B_vectors2 is not None:
        B_components2 = B_vectors2.transpose()
        ax1.plot(
            times,
            B_components2[0],
            color=color2,
            label=label2,
            linestyle=linestyle2,
            linewidth=linewidth2,
            alpha=.7
        )
        ax2.plot(
            times,
            B_components2[1],
            color=color2,
            linestyle=linestyle2,
            linewidth=linewidth2,
            alpha=.7

        )
    else:
        pass
    if B_vectors3 is not None:
        B_components3 = B_vectors3.transpose()
        ax1.plot(
            times,
            B_components3[0],
            color=color3,
            label=label3,
            linestyle=linestyle3,
            linewidth=linewidth2,
            alpha=alpha
        )
        ax2.plot(
            times,
            B_components3[1],
            color=color3,
            linestyle=linestyle3,
            linewidth=linewidth2,
            alpha=alpha
        )
    else:
        pass
    fig.text(0.53, 0.01, "Time (hours)", ha="center", fontsize=20)
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(
        lines,
        labels,
        loc="upper center",
        ncols=3,
        bbox_to_anchor=bbox_coords,
        frameon=False,
    )
    fig.tight_layout(pad=1.5, w_pad=0.5, h_pad=1.5)

    if filename is not None:
        plt.savefig(
            f"/Users/mikoo/OneDrive - Imperial College London/Year 4/project/plots/{filename}.png"
        )
    else:
        pass

    return plt.show()


def plot_B_evolution_box(
    B_vectors1,
    times,
    B_vectors2: np.ndarray = None,
    B_vectors3: np.ndarray = None,
    label1: str = None,
    label2: str = None,
    label3: str = None,
    filename: str = None,
    color1: str = None,
    color2: str = None,
    linestyle1: str = None,
    linestyle2: str = None,
    linestyle3: str = None,
    linewidth1: float = None,
    linewidth2: float = None,
    bbox_coords: tuple = None,
    pad: float = None,
):
    B_components = B_vectors1.transpose()
    times = times / 3600  # convert to hours
    params = {
        "axes.labelsize": 20,
        "font.size": 18,
        "legend.fontsize": 22,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "figure.figsize": [14, 8],
    }
    plt.rcParams.update(params)

    Bx = B_components[0]
    By = B_components[1]
    Bz = B_components[2]
    Bmag = np.linalg.norm(B_components, axis=0)

    fig, axs = plt.subplots(2, 2, sharex=True)

    ax1 = axs[0, 0]
    ax1.plot(
        times,
        Bx,
        label=label1,
        color=color1,
        linewidth=linewidth1,
        linestyle=linestyle1,
    )
    ax1.set_ylabel("Bx (nT)")
    ax1.set_xlim(times[0], times[-1])
    ax1.grid()
    ax1.minorticks_on()

    ax2 = axs[0, 1]
    ax2.plot(times, By, color=color1, linewidth=linewidth1, linestyle=linestyle1)
    ax2.set_ylabel("By (nT)")
    ax2.set_xlim(times[0], times[-1])
    ax2.grid()
    ax2.minorticks_on()

    ax3 = axs[1, 0]
    ax3.plot(times, Bz, color=color1, linewidth=linewidth1, linestyle=linestyle1)
    ax3.set_ylabel("Bz (nT)")
    ax3.set_xlim(times[0], times[-1])
    # ax3.set_xlabel("Time (hours)")
    ax3.set_xlabel("Time (seconds)")

    ax3.grid()
    ax3.minorticks_on()

    ax4 = axs[1, 1]
    ax4.plot(times, Bmag, color=color1, linewidth=linewidth1, linestyle=linestyle1)
    ax4.set_ylabel(r"$\vert$B$\vert$ (nT)")
    ax4.set_xlabel("Time (hours)")
    ax4.set_xlabel("Time (seconds)")
    ax4.set_xlim(times[0], times[-1])
    ax4.grid()
    ax4.minorticks_on()

    if B_vectors2 is not None:
        B_components2 = B_vectors2.transpose()
        ax1.plot(
            times,
            B_components2[0],
            color=color2,
            label=label2,
            linestyle=linestyle2,
            linewidth=linewidth2,
        )
        ax2.plot(
            times,
            B_components2[1],
            color=color2,
            linestyle=linestyle2,
            linewidth=linewidth2,
        )
        ax3.plot(
            times,
            B_components2[2],
            color=color2,
            linestyle=linestyle2,
            linewidth=linewidth2,
        )
        ax4.plot(
            times,
            np.linalg.norm(B_components2, axis=0),
            color=color2,
            linestyle=linestyle2,
            linewidth=linewidth2,
        )
    else:
        pass
    if B_vectors3 is not None:
        B_components3 = B_vectors3.transpose()
        ax1.plot(
            times,
            B_components3[0],
            color=color2,
            label=label3,
            linestyle=linestyle3,
            linewidth=linewidth2,
        )
        ax2.plot(
            times,
            B_components3[1],
            color=color2,
            linestyle=linestyle3,
            linewidth=linewidth2,
        )
        ax3.plot(
            times,
            B_components3[2],
            color=color2,
            linestyle=linestyle3,
            linewidth=linewidth2,
        )
        ax4.plot(
            times,
            np.linalg.norm(B_components3, axis=0),
            color=color2,
            linestyle=linestyle3,
            linewidth=linewidth2,
        )
    else:
        pass

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(
        lines,
        labels,
        loc="upper center",
        ncols=3,
        bbox_to_anchor=bbox_coords,
        frameon=False,
    )
    fig.tight_layout(pad=pad, w_pad=1, h_pad=0.5)

    if filename is not None:
        plt.savefig(
            f"/Users/mikoo/OneDrive - Imperial College London/Year 4/project/plots/{filename}.png"
        )
    else:
        pass

    return plt.show()


def plot_B_evolution_column_3(
    z,
    rho,
    B_vectors1,
    times,
    B_vectors2: np.ndarray = None,
    B_vectors3: np.ndarray = None,
    label1: str = None,
    label2: str = None,
    label3: str = None,
    filename: str = None,
    color1: str = None,
    color2: str = None,
    linestyle1: str = None,
    linestyle2: str = None,
    linestyle3: str = None,
    linewidth1: float = None,
    linewidth2: float = None,
    bbox_coords: tuple = None,
    pad: float = None,
):
    B_components = B_vectors1.transpose()
    times = times / 3600  # convert to hours
    params = {
        "axes.labelsize": 20,
        "font.size": 18,
        "legend.fontsize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "figure.figsize": [14, 10],
    }
    plt.rcParams.update(params)

    By = B_components[1]
    Bz = B_components[2]
    Bmag = np.linalg.norm(B_components, axis=0)

    fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True)

    ax0 = axs[0]
    ax0.plot(
        times,
        rho,
        label=r"$\rho$, z",
        color=color1,
        linewidth=linewidth1,
        linestyle="-.",
    )
    ax0.set_ylabel(r"$\rho$ (RJ)")
    ax0.set_xlim(times[0], times[-1])
    ax0.grid()
    ax0.minorticks_on()

    ax1 = axs[1]
    ax1.plot(times, z, color=color1, linewidth=linewidth1, linestyle="-.")
    ax1.set_ylabel(r"z (RJ)")
    ax1.set_xlim(times[0], times[-1])
    ax1.grid()
    ax1.minorticks_on()

    ax2 = axs[2]
    ax3 = axs[3]
    ax4 = axs[4]

    if B_vectors2 is not None:
        B_components2 = B_vectors2.transpose()
      # ax1.plot(times, B_components2[0], label=label2, color=color2, linestyle=linestyle2, linewidth=linewidth2)
        ax2.plot(
        times,
        B_components2[1],
        color=color2,
        label=label2,
        linestyle=linestyle2,
        linewidth=linewidth2,
        )
        ax3.plot(
            times,
            B_components2[2],
            color=color2,
            linestyle=linestyle2,
            linewidth=linewidth2,
        )
        ax4.plot(
            times,
            np.linalg.norm(B_components2, axis=0),
            color=color2,
            linestyle=linestyle2,
            linewidth=linewidth2,
        )
    else:
        pass

    if B_vectors3 is not None:
        B_components3 = B_vectors3.transpose()
        # ax1.plot(times, B_components3[0], label=label3, color=color2, linestyle=linestyle3, linewidth=linewidth2)
        ax2.plot(
            times,
            B_components3[1],
            color=color2,
            label=label3,
            linestyle=linestyle3,
            linewidth=linewidth2,
        )
        ax3.plot(
            times,
            B_components3[2],
            color=color2,
            linestyle=linestyle3,
            linewidth=linewidth2,
        )
        ax4.plot(
            times,
            np.linalg.norm(B_components3, axis=0),
            color=color2,
            linestyle=linestyle3,
            linewidth=linewidth2,
        )

        ax2.plot(
            times,
            By,
            label=label1,
            color=color1,
            linewidth=linewidth1,
            linestyle=linestyle1,
        )
        ax2.set_ylabel(r"By (nT)")
        ax2.set_xlim(times[0], times[-1])
        ax2.grid()
        ax2.minorticks_on()

        ax3.plot(times, Bz, color=color1, linewidth=linewidth1, linestyle=linestyle1)
        ax3.set_ylabel("Bz (nT)")
        ax3.set_xlim(times[0], times[-1])
        ax3.grid()
        ax3.minorticks_on()

        ax4.plot(times, Bmag, color=color1, linewidth=linewidth1, linestyle=linestyle1)
        ax4.set_ylabel(r"$\vert$B$\vert$ (nT)")
        ax4.set_xlim(times[0], times[-1])
        ax4.set_xlabel("Time (hours)")
        ax4.grid()
        ax4.minorticks_on()

        ax2.fill_between(
            x=times, y1=B_components2[1], y2=B_components[1], alpha=0.3, color="grey"
        )
        ax3.fill_between(
            x=times, y1=B_components2[2], y2=B_components[2], alpha=0.3, color="grey"
        )
        ax4.fill_between(
            x=times,
            y1=np.linalg.norm(B_components2, axis=0),
            y2=np.linalg.norm(B_components, axis=0),
            alpha=0.3,
            color="grey",
        )
    else:
        pass

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(
        lines,
        labels,
        loc="upper center",
        ncols=4,
        bbox_to_anchor=bbox_coords,
        frameon=False,
    )
    fig.tight_layout(pad=pad, w_pad=1, h_pad=0.5)

    if filename is not None:
        plt.savefig(
            f"/Users/mikoo/OneDrive - Imperial College London/Year 4/project/plots/{filename}.png"
        )
    else:
        pass

    return plt.show()



