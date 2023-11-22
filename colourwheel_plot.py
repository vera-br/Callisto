from functions import *
from adjustText import adjust_text

# ----------- GALILEO ------------

galileo_wrt_callisto_cphio_CA = closest_approach_data_G("galileo", "callisto", "cphio", "G")
callisto_wrt_jupiter_cphio_CA = closest_approach_data_G("callisto", "jupiter", "cphio", "G")
G_callisto_wrt_jupiter_JSO_CA = closest_approach_data_G("callisto", "jupiter", "jupsunorb", "G")
sun_wrt_callisto_cphio_CA = closest_approach_data_G("sun", "callisto", "cphio", "G")
callisto_wrt_jupiter_SIII_mag_CA = closest_approach_data_G("callisto", "jupiter", "SIII_mag", "G")
jupiter_wrt_sun_IAU_CA = closest_approach_data_G('jupiter', 'sun', 'IAU_SUN', 'G')

# get jupiter-sun angles
G_zenith = []
for orbit, vector in jupiter_wrt_sun_IAU_CA.items():
    G_zenith.append(vector[6])

# get sun-juice angle
G_sun_cphio_vector = []
for orbit, vector in sun_wrt_callisto_cphio_CA.items():
    G_sun_cphio_vector.append(vector[1:4])

galileo_cphio_vector = []
for orbit, vector in galileo_wrt_callisto_cphio_CA.items():
    galileo_cphio_vector.append(vector[1:4])

sun_galileo_angle = []
for i in range(len(G_sun_cphio_vector)):
    angle = angle_between(G_sun_cphio_vector[i], galileo_cphio_vector[i])
    sun_galileo_angle.append(np.degrees(angle))

# ----------- JUICE ------------

juice_wrt_callisto_cphio_CA = get_closest_approach_data("juice", "callisto", "cphio", "J")
callisto_wrt_jupiter_cphio_CA = get_closest_approach_data("callisto", "jupiter", "cphio", "J")
J_callisto_wrt_jupiter_JSO_CA = get_closest_approach_data("callisto", "jupiter", "jupsunorb", "J")
sun_wrt_callisto_cphio_CA = get_closest_approach_data("sun", "callisto", "cphio", "J")
callisto_wrt_jupiter_SIII_mag_CA = get_closest_approach_data("callisto", "jupiter", "SIII_mag", "J")
jupiter_wrt_sun_IAU_CA = get_closest_approach_data('jupiter', 'sun', 'IAU_SUN', 'J')

# get jupiter-sun angles
J_zenith = []
for orbit, vector in jupiter_wrt_sun_IAU_CA.items():
    J_zenith.append(vector[6])

# get sun-juice angle
J_sun_cphio_vector = []
for orbit, vector in sun_wrt_callisto_cphio_CA.items():
    J_sun_cphio_vector.append(vector[1:4])

juice_cphio_vector = []
for orbit, vector in juice_wrt_callisto_cphio_CA.items():
    juice_cphio_vector.append(vector[1:4])

sun_juice_angle = []
for i in range(len(J_sun_cphio_vector)):
    angle = angle_between(J_sun_cphio_vector[i], juice_cphio_vector[i])
    sun_juice_angle.append(np.degrees(angle))


# Data for the colour wheel
xval = np.arange(-np.pi, np.pi, 0.01)
yval = np.ones_like(xval)

colormap = plt.get_cmap('hsv_r')
norm = mpl.colors.Normalize(-np.pi, np.pi)


# Plot
fig = plt.figure(figsize=(10,5))
ax1 = plt.subplot(121)

# Defining and applying the common limits
lim = 30

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
for orbit, vector in J_callisto_wrt_jupiter_JSO_CA.items():

    if abs(sun_juice_angle[i]) < 90:
        s = ax1.scatter(vector[1] / R_J, vector[2] / R_J, c=J_zenith[i], vmin=-np.pi, vmax=np.pi, s=80, cmap=colormap, marker='*')
    else:
        s = ax1.scatter(vector[1] / R_J, vector[2] / R_J, c=J_zenith[i], vmin=-np.pi, vmax=np.pi, s=30, cmap=colormap)

    i += 1
    texts.append(ax1.text(vector[1] / R_J, vector[2] / R_J, 'J%s' % i)) # save orbit labels in list

i=0
for orbit, vector in G_callisto_wrt_jupiter_JSO_CA.items():

    if abs(sun_galileo_angle[i]) < 90:
        s = ax1.scatter(vector[1] / R_J, vector[2] / R_J, c=G_zenith[i], vmin=-np.pi, vmax=np.pi, cmap=colormap, s=80, marker='*', facecolor='none')
    else:
        s = ax1.scatter(vector[1] / R_J, vector[2] / R_J, c=G_zenith[i], vmin=-np.pi, vmax=np.pi, cmap=colormap, s=30, marker='o')

    i += 1
    texts.append(ax1.text(vector[1] / R_J, vector[2] / R_J, 'G%s' % i)) # save orbit labels in list

# Adjust the positions of annotations so they don't overlap
adjust_text(texts)

# colourbar
#cbar = fig.colorbar(s)
#cbar.set_title('Jupiter-Sun Angle in IAU_SUN [degrees]')

legend = ax1.legend(handles=[dayside_marker, nightside_marker], loc='lower left')
ax1.set_title('Closest Approaches in JSO')

# ------------ COLOUR WHEEL --------------

ax2 = plt.subplot(122, projection = 'polar')
ax2.scatter(xval, yval, c=xval, s=300, cmap=colormap, norm=norm, linewidths=0)
ax2.scatter(0,0, color='gold')

labels =[]
for i in range(len(J_zenith)):
    if sun_juice_angle[i] < 90:
        ax2.scatter(J_zenith[i], 1, color='black', marker='*')
    else:
        ax2.scatter(J_zenith[i], 1, color='black', s=10)
    labels.append(ax2.text(J_zenith[i], 1, 'J%s' % (i+1)))

for i in range(len(G_zenith)):
    if sun_galileo_angle[i] < 90:
        ax2.scatter(G_zenith[i], 1, color='black', marker='*', facecolors='none')
    else:
        ax2.scatter(G_zenith[i], 1, color='black', marker='o')
    labels.append(ax2.text(G_zenith[i], 1, 'G%s' % (i+1)))

adjust_text(labels)
ax2.set_yticks([])
ax2.set_title('Jupiter-Sun Angle in IAU_SUN')

plt.tight_layout()
plt.show()