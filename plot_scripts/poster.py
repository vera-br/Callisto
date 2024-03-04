# plots for inducing field on poster
from jupiter_field import *
from khurana1997 import *
from plot_scripts.plot_field_components import *

# load data
callisto_jupiter_SIII = get_spice_data('callisto', 'jupiter', 'SIII', 'J')
callisto_jupiter_JSO = get_spice_data('callisto', 'jupiter', 'jupsunorb', 'J')
callisto_jupiter_SIII_mag = get_spice_data('callisto', 'jupiter', 'SIII_mag', 'J')

flyby_n = 2

orbit_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]

# create array with time series and constant position in CphiO
time = np.linspace(0, 24*3600, len(orbit_SIII[0]))

constant_pos = np.array([0, 150e3, 0]) #200km above surface on y axis
constant_pos = constant_pos.reshape(1, -1) 
pos_series = np.repeat(constant_pos, len(time), axis=0)

orbit_200km = np.c_[time, pos_series]


B_Jupiter = Bext_Community(orbit_SIII)
B_sheet = B_sheet_khurana2(orbit_JSO, orbit_SIII_mag)

#-------------

B_field_mag_jup = (B_Jupiter[:,0]**2 + B_Jupiter[:,1]**2 + B_Jupiter[:,2]**2)**0.5
B_field_mag_jup = pd.Series(dict(zip(time, B_field_mag_jup)))

B_field_mag_sheet = (B_sheet[:,0]**2 + B_sheet[:,1]**2 + B_sheet[:,2]**2)**0.5
B_field_mag_sheet = pd.Series(dict(zip(time, B_field_mag_sheet)))


fig = plt.figure(figsize=(8,4), dpi=300, constrained_layout=True)
ax = fig.gca()

plt.plot(time/3600, B_field_mag_sheet, label="Magnetodisk Field", color="#dc267f")
plt.plot(time/3600, B_field_mag_jup, label="Jupiter's Field", color="#648fff")

ax.tick_params(axis='both', direction='in',top = True, right = True, which='both', labelsize=14)
ax.set_xticks(np.arange(0, 29, step=6))
#ax.set_yticks()
ax.yaxis.set_minor_locator(AutoMinorLocator()) 
ax.xaxis.set_minor_locator(AutoMinorLocator())

ax.set_xlim(0,24)

plt.legend(loc='upper right', fontsize=14)
#ax.grid(color='xkcd:dark blue',alpha =0.2)

ax.set_xlabel("Time [hours]", fontsize=16)
ax.set_ylabel("|B| [nT]", fontsize=16)
#ax.set_title(field_type + " field %skm above Callisto's surface" % (altitude))
#plt.tight_layout()
plt.savefig("variation.png", facecolor=('b', 0), pad_inches=0.1)
