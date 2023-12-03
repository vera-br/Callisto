"""
The following script calculates the magnetic field components for a given trajectory
"""
from jupiter_external import *
from induced_field import *
from plot_scripts.plot_components import *
from plot_scripts.spatial_variations import position_evolution_column
from current_sheet import B_sheet_mag



# constants
pi = constants.pi
mu0 = constants.mu_0
RJ = 71492e3  # Jupiter radius
RC = 2410.1e3  # callisto radius

### Load in data from SPICE kernels, data loaded in has units km and km/s
# data paths
juice_callisto_data_path = "./spice_data/juice_wrt_callisto_cphio_J4.csv"
juice_jupiter_data_path = "./spice_data/juice_wrt_jupiter_cphio_J4.csv"
callisto_jupiter_data_path = "./spice_data/callisto_wrt_jupiter_cphio_J4.csv"
juice_jupiter_mag_data_path = "./spice_data/juice_wrt_jupiter_SIII_mag_J4.csv"

# juice to callisto
juice_wrt_callisto = np.loadtxt(
    juice_callisto_data_path, delimiter=",", unpack=True
)
juice_wrt_callisto = juice_wrt_callisto.transpose()
Xjg, Yjg, Zjg, vxjg, vyjg, vzjg, t = juice_wrt_callisto

# juice to jupiter
juice_wrt_jupiter = np.loadtxt(
    juice_jupiter_data_path, delimiter=",", unpack=True
)
juice_wrt_jupiter = juice_wrt_jupiter.transpose()
Xjj, Yjj, Zjj, vxjj, vyjj, vzjj, _ = juice_wrt_jupiter

# callisto to jupiter
callisto_wrt_jupiter = np.loadtxt(
    callisto_jupiter_data_path, delimiter=",", unpack=True
)
callisto_wrt_jupiter = callisto_wrt_jupiter.transpose()
Xgj, Ygj, Zgj, vxgj, vygj, vzgj, _ = callisto_wrt_jupiter

# 
juice_wrt_jupiter_mag = np.loadtxt(
    juice_jupiter_mag_data_path, delimiter=",", unpack=True
)
juice_wrt_jupiter_mag = juice_wrt_jupiter_mag.transpose()
Xjjm, Yjjm, Zjjm, vxjjm, vyjjm, vzjjm, _ = juice_wrt_jupiter_mag


# convert positions to m
juice_wrt_callisto = 1e3 * np.array([Xjg, Yjg, Zjg]).transpose()
juice_wrt_jupiter = 1e3 * np.array([Xjj, Yjj, Zjj]).transpose()
callisto_wrt_jupiter = 1e3 * np.array([Xgj, Ygj, Zgj]).transpose()

# set initial time to zero
t = t - t[0]

# induced field parameters
rm = 2410.1e3  # callisto radius
r0 = rm - 100e3  # outer conducting layer radius
r1 = rm - 200e3  # inner conducting layer radius

# rotation speeds
J_spin_period = 10.1 * 3600  # actual period is 9.9 hours
C_orbit_period = 17 * 24 * 3600
C_omega = 2 * pi / C_orbit_period  # hours^-1
J_omega = 2 * pi / J_spin_period  # hours^-1

# induced bessel parameters
omega_synodic = 2 * pi / J_spin_period
omega_eccentric = 2 * pi / C_orbit_period

# parameters for cylindrical plasma sheet
R0 = 7.8  # disc inner radius (RJ)
R1 = 51.4  # disc outer radius (RJ)
D = 3.6  # disc half thickness (RJ)
Icon = 139.6  # current constant = mu0 * I / 2 (nT)
thetaD = np.radians(9.3)  # disc normal from rotation axis (radians)
phiD = np.radians(204.2)  # azimuth angle of disc normal (radians)

# calculate magnetic fields
# B_intrinsic = B_intrinsic_vectors_dipole(positions=juice_wrt_ganymede)
B_external = Bext_full(positions=juice_wrt_jupiter, times=t)

Bext_ecc = Bext_eccentric_variation(positions=callisto_wrt_jupiter)
Bind_ecc_superconductor = B_induced_superconductor(
    pos_vectors=juice_wrt_callisto, Bext_vectors=Bext_ecc, r0=r0, rm=rm
)
Bind_ecc_finite = B_induced_finite_conductivity(
    pos_vectors=juice_wrt_callisto,
    Bext_vectors=Bext_ecc,
    sigma=0.03,
    omega=omega_synodic,
    rm=rm,
    r0=r0,
    r1=rm - 150e3 - 300e3,
)

Bext_syn = Bext_synodic_variation(pos_avg=np.mean(juice_wrt_jupiter, axis=0), times=t)

Bind_syn_superconductor = B_induced_superconductor(
    pos_vectors=juice_wrt_callisto, Bext_vectors=Bext_syn, r0=rm - 150e3, rm=rm
)
Bind_syn_finite = B_induced_finite_conductivity(
    pos_vectors=juice_wrt_callisto,
    Bext_vectors=Bext_syn,
    sigma=0.03,
    omega=omega_synodic,
    rm=rm,
    r0=r0,
    r1=rm - 150e3 - 300e3,
)

pos = [Xjjm / RJ, Yjjm / RJ, Zjjm /RJ]
rho, z, Bsheet = B_sheet_mag(
    positions=pos, R0=R0, R1=R1, D=D, I_constant=Icon
)


###
# plot sheet fields for varying D
# plot_B_evolution_column_3(
#     B_vectors1=B_plasma_sheet3,
#     B_vectors2=B_plasma_sheet2,
#     B_vectors3=B_plasma_sheet1,
#     z=z, rho=rho,
#     times=t,
#     color1="black",
#     color2="black",
#     linestyle1="--",
#     linestyle2="dotted",
#     linestyle3="-",
#     linewidth1=1.,
#     linewidth2=1.,
#     label2=r"D = 2.6 RJ",
#     label1=r"D = 4.6 RJ",
#     label3=r"D = 3.6 RJ",
#     bbox_coords=(0.53, 1.005),
#     pad=2,
#     filename="gco_500/B_sheet_comps"
# )

# juice_trajectory(positions=juice_wrt_callisto.transpose()/RG)
# plot_B_evolution_column(B_vectors1=B_external, times=t, filename="Bext_comps", color1="black", linewidth1=1.)

# position_evolution_column(B_vectors1=juice_wrt_callisto / RG, times=t, color="black")

###
# Plot ECCENTRIC FIELD VARIATION
# plot_B_evolution_column(
#     B_vectors1=Bext_ecc_var,
#     # B_vectors2=Bind_ecc_superconductor,
#     B_vectors2=Bind_ecc_finite,
#     times=t,
#     filename="Bind_eccentric_finite",
#     color1="black",
#     color2="black",
#     linestyle1="dotted",
#     linestyle2="-",
#     linestyle3="--",
#     linewidth1=1.5,
#     linewidth2=0.9,
#     label1=r"$\delta$B",
#     # label2=r"B$_{\mathrm{ind}}$($\sigma\rightarrow\infty$)",
#     label2=r"B$_{\mathrm{ind}}$($\sigma=$2 S/m)",
#     bbox_coords=(0.53, 1.01),
#     pad=2
# )
###

###
# plot SYNODIC FIELD VARIATION POSTER
# plot_B_evolution_poster_ind(
#     B_vectors1=Bext_syn_var,
#     B_vectors2=Bind_syn_superconductor,
#     B_vectors3=Bind_syn_finite,
#     times=t,
#     filename="Bind_sigma_poster",
#     color1="black",
#     color2="blue",
#     color3="magenta",
#     linestyle1="-",
#     linestyle2="-",
#     linestyle3="-",
#     linewidth1=2,
#     linewidth2=2,
#     label1=r"$\mathrm{\delta}$B$_{\mathrm{ext}}$",
#     label2=r"B$_{\mathrm{ind}}$($\mathrm{\sigma\rightarrow\infty}$)",
#     label3=r"B$_{\mathrm{ind}}$($\mathrm{\sigma=}$0.03 S/m)",
#     # label2=r"B$_{\mathrm{ind}}$($h$ = 300 km)",
#     # label3=r"B$_{\mathrm{ind}}$($h$ = 50 km)",
#     bbox_coords=(0.535, 1.042),
#     pad=1.3
# ) 
###

###
# plot TOTAL FIELD AND POWER

Btot_super = B_external + Bind_syn_superconductor + Bind_ecc_superconductor #+ B_plasma_sheet1 + B_intrinsic
Btot_finite = B_external + Bind_syn_finite + Bind_ecc_finite #+ B_plasma_sheet1 + B_intrinsic
plot_B_evolution_column_3(
    z=z,
    rho=rho,
    B_vectors1=Btot_finite,
    times=t,
    filename="Btotal_super",
    # filename="Btotal_finite",
    color1="black",
    color2="black",
    linestyle1="-",
    linewidth1=.9,
    bbox_coords=(0.53, 1.01),
    pad=2
)
# fourier_decomposition(B_vectors=Btot_finite, time=t, peak_sens=0.00001, color="black", filename="gco_500/Btot_fourier")
###

###
# Bz induction in eccentric orbit plot
# params = {
#         "axes.labelsize": 14,
#         "font.size": 18,
#         "legend.fontsize": 22,
#         "xtick.labelsize": 11,
#         "ytick.labelsize": 11,
#         "figure.figsize": [6, 6],
#     }
# plt.rcParams.update(params)
# fig, axs = plt.subplots(1, 1, sharex=True)
# axs.plot(t / 3600, Bind_ecc_superconductor.transpose()[-1], color="black", linestyle="-", linewidth=1.)
# axs.plot(t / 3600, Bext_ecc_var.transpose()[-1], color="black", linestyle="-.", linewidth=1.)
# axs.set_xlim(0, 172.57)
# axs.set_xlabel("Time (hours)")
# axs.set_ylabel("Bz (nT)")
# axs.grid()
# axs.minorticks_on()
# fig.tight_layout(pad=1, w_pad=0.5, h_pad=.5)
# plt.savefig(f"/Users/ciaranjones/Documents/University/Year4/MSci_Project/Plots/report_plots/gco_500/ecc_poster.pdf")
# plt.show()
###

###
# ECC + SYN POSTER
# params = {
#     "axes.labelsize": 20,
#     "font.size": 18,
#     "legend.fontsize": 20,
#     "xtick.labelsize": 18,
#     "ytick.labelsize": 18,
#     "figure.figsize": [12, 6],
# }
# plt.rcParams.update(params)
# fig, axs = plt.subplots(nrows=1, ncols=2)
# ax1 = axs[0]
# ax2 = axs[1]
#
# ax1.plot(t / 3600, Bext_ecc_var.transpose()[-1], color="black", linestyle="-", linewidth=2)
#
# ax1.plot(t / 3600, Bind_ecc_superconductor.transpose()[-1], color="blue", linestyle="-", linewidth=2, alpha=0.5)
# ax1.set_xlim(0, 172.57)
# ax1.set_ylim(-.8, .8)
# ax1.set_ylabel("Bz (nT)")
# ax1.grid()
# ax1.minorticks_on()
#
#
# ax2.plot(t / 3600, np.linalg.norm(Bext_syn_var.transpose(), axis=0), color="black", linestyle="-", linewidth=2)
# ax2.plot(t / 3600, np.linalg.norm(Bind_ecc_superconductor.transpose() + Bind_syn_superconductor.transpose(), axis=0), color="blue", alpha=0.7, linestyle="-", linewidth=2)
# ax2.plot(t / 3600, np.linalg.norm(Bind_ecc_finite.transpose() + Bind_syn_finite.transpose(), axis=0), color="magenta", alpha=.7, linestyle="-", linewidth=2)
# ax2.set_xlim(0, 2 * 10.53)
# ax2.set_ylim(0, 30)
# ax2.set_ylabel(r"$\vert$B$\vert$ (nT)")
# ax2.grid()
# ax2.minorticks_on()
#
# fig.text(0.53, 0.01, 'Time (hours)', ha='center', fontsize=20)
# fig.tight_layout(pad=1.2, w_pad=1.5, h_pad=.5)
# plt.savefig(f"/Users/ciaranjones/Documents/University/Year4/MSci_Project/Plots/report_plots/gco_500/ecc_syn_poster.pdf")
# plt.show()
