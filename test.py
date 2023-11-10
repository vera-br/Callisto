from functions import *
from itertools import cycle
from adjustText import adjust_text


# load orbit data
juice_wrt_callisto_cphio = get_spice_data('juice', 'callisto', 'cphio')
callisto_wrt_jupiter_cphio = get_spice_data("callisto", "jupiter", "jupsunorb")
sun_wrt_callisto_cphio = get_spice_data("sun", "callisto", "cphio")
callisto_wrt_jupiter_SIII = get_spice_data("callisto", "jupiter", "SIII")
jupiter_wrt_sun_IAU = get_spice_data('jupiter', 'sun', 'IAU_SUN')

Galileo, Galileo_meas = get_pds_data()

# save CA data
juice_wrt_callisto_cphio_CA, sun_wrt_callisto_cphio_CA, callisto_wrt_jupiter_SIII_CA, jupiter_wrt_sun_IAU_CA = closest_approach_data_4(juice_wrt_callisto_cphio, sun_wrt_callisto_cphio, callisto_wrt_jupiter_SIII, jupiter_wrt_sun_IAU)
Galileo_CA = closest_approach_data(Galileo)


# 3D Plot
fig = plt.figure()
ax = fig.add_subplot()

# Defining and applying the common limits
lim = 35

common_xlim = (-lim, lim)  
common_ylim = (-lim, lim) 

ax.set_xlim(common_xlim)
ax.set_ylim(common_ylim)

# Set axis labels
ax.set_xlabel('x [$R_J$]')
ax.set_ylabel('y [$R_J$]')

# plot jupiter
jup = plt.Circle((0, 0), 1, color='black')
ax.add_patch(jup)

# get jupiter-sun angles
azimuthal = []
for orbit, vector in jupiter_wrt_sun_IAU_CA.items():
    azimuthal.append(np.degrees(vector[6]))

# get sun-juice angle
sun_cphio_vector = []
for orbit, vector in sun_wrt_callisto_cphio_CA.items():
    sun_cphio_angle.append(vector[1:4])

juice_cphio_vector = []
for orbit, vector in juice_wrt_callisto_cphio.items():
    juice_cphio_angle.append(vector[1:4])

juice_sun_angle = angle_between(juice_cphio_vector, sun_cphio_vector)

texts = []
i = 0

for orbit, vector in callisto_wrt_jupiter_SIII_CA.items():

    if abs(sun-cal-angle[i]) > 90:
        s = ax.scatter(vector[1] / R_J, vector[2] / R_J, c=azimuthal[i], vmin=min(azimuthal), vmax=max(azimuthal), s=20, cmap='plasma', marker='*')
    else
        s = ax.scatter(vector[1] / R_J, vector[2] / R_J, c=azimuthal[i], vmin=min(azimuthal), vmax=max(azimuthal), s=20, cmap='plasma', marker='.')

    i += 1
    texts.append(ax.text(vector[1] / R_J, vector[2] / R_J, '%s' % i)) # save orbit labels in list

# Adjust the positions of annotations
adjust_text(texts)

cbar = fig.colorbar(s)
cbar.set_label(' Jupiter-Sun Angle in IAU_SUN [degrees]')
plt.title('Closest Approaches in SIII')
plt.tight_layout()
plt.show()