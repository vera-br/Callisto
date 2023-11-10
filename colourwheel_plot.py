from functions import *
from adjustText import adjust_text


# load orbit data
juice_wrt_callisto_cphio = get_spice_data('juice', 'callisto', 'cphio')
callisto_wrt_jupiter_cphio = get_spice_data("callisto", "jupiter", "cphio")
callisto_wrt_jupiter_JSO = get_spice_data("callisto", "jupiter", "jupsunorb")
sun_wrt_callisto_cphio = get_spice_data("sun", "callisto", "cphio")
callisto_wrt_jupiter_SIII = get_spice_data("callisto", "jupiter", "SIII")
jupiter_wrt_sun_IAU = get_spice_data('jupiter', 'sun', 'IAU_SUN')

Galileo, Galileo_meas = get_pds_data()

# save CA data
juice_wrt_callisto_cphio_CA, sun_wrt_callisto_cphio_CA, callisto_wrt_jupiter_JSO_CA, jupiter_wrt_sun_IAU_CA = closest_approach_data_4(juice_wrt_callisto_cphio, sun_wrt_callisto_cphio, callisto_wrt_jupiter_JSO, jupiter_wrt_sun_IAU)
Galileo_CA = closest_approach_data(Galileo)


# get jupiter-sun angles
azimuthal = []
for orbit, vector in jupiter_wrt_sun_IAU_CA.items():
    azimuthal.append(np.degrees(vector[6]))

# get sun-juice angle
sun_cphio_vector = []
for orbit, vector in sun_wrt_callisto_cphio_CA.items():
    sun_cphio_vector.append(vector[1:4])

juice_cphio_vector = []
for orbit, vector in juice_wrt_callisto_cphio_CA.items():
    juice_cphio_vector.append(vector[1:4])

sun_juice_angle = []
for i in range(len(sun_cphio_vector)):
    angle = angle_between_vectors(sun_cphio_vector[i], juice_cphio_vector[i])
    sun_juice_angle.append(angle)


# Data for the colour wheel
xval = np.arange(-np.pi, np.pi, 0.01)
yval = np.ones_like(xval)

colormap = plt.get_cmap('hsv_r')
norm = mpl.colors.Normalize(-np.pi, np.pi)


# Plot
fig = plt.figure(figsize=(12,6))
ax1 = plt.subplot(121)

# Defining and applying the common limits
lim = 28

common_xlim = (-lim, lim)  
common_ylim = (-lim, lim) 

ax1.set_xlim(common_xlim)
ax1.set_ylim(common_ylim)

# Set axis labels
ax1.set_xlabel('x [$R_J$]')
ax1.set_ylabel('y [$R_J$]')

# plot jupiter
jup = plt.Circle((0, 0), 1, color='xkcd:dull brown')
ax1.add_patch(jup)
ax1.annotate('Jupiter', (-2,-3))

dayside_marker = mlines.Line2D([], [], color='black', marker='*', linestyle='None', markersize=10, label='Dayside')
nightside_marker = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=7, label='Nightside')

texts = []
i = 0

for orbit, vector in callisto_wrt_jupiter_JSO_CA.items():

    if sun_juice_angle[i] > 90:
        s = ax1.scatter(vector[1] / R_J, vector[2] / R_J, c=azimuthal[i], vmin=min(azimuthal), vmax=max(azimuthal), s=80, cmap=colormap, marker='*')
    else:
        s = ax1.scatter(vector[1] / R_J, vector[2] / R_J, c=azimuthal[i], vmin=min(azimuthal), vmax=max(azimuthal), s=40, cmap=colormap)

    i += 1
    texts.append(ax1.text(vector[1] / R_J, vector[2] / R_J, 'C%s' % i)) # save orbit labels in list

# Adjust the positions of annotations so they don't overlap
adjust_text(texts)

# colourbar
#cbar = fig.colorbar(s)
#cbar.set_title('Jupiter-Sun Angle in IAU_SUN [degrees]')

legend = ax1.legend(handles=[dayside_marker, nightside_marker], loc='lower left')
ax1.set_title('Closest Approaches in JSO')

ax2 = plt.subplot(122, projection = 'polar')
ax2.scatter(xval, yval, c=xval, s=300, cmap=colormap, norm=norm, linewidths=0)
ax2.set_yticks([])
ax2.set_title('Jupiter-Sun Angle in IAU_SUN')

plt.tight_layout()
plt.show()