# plot script for field components evolving in time

import pandas as pd
from pandas import Timestamp
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

# define constants
R_C = 2410.3 * 1e3
R_J = 71492 * 1e3
R_S = 696340 * 1e3
AU = 150000000 * 1e3


def plot_time_evolution(B_field, orbit_cphio, orbit_CA, flyby_n, field_type):

    # convert time into Timestamp format
    OFFSET = datetime(2000,1,1,12) - datetime(1970,1,1) # difference between J2000 and UTC

    time = orbit_cphio[0]
    time = [Timestamp((datetime.utcfromtimestamp(timestamp) + OFFSET).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time]
    time_CA = Timestamp((datetime.utcfromtimestamp(orbit_CA[0]) + OFFSET).strftime('%Y-%m-%d %H:%M:%S'))

    # radial distance of juice wrt callisto
    radial = pd.Series(dict(zip(time, orbit_cphio[4] / R_C)))

    B_field_x = pd.Series(dict(zip(time, B_field[:,0])))
    B_field_y = pd.Series(dict(zip(time, B_field[:,1])))
    B_field_z = pd.Series(dict(zip(time, B_field[:,2])))

    B_field_mag = (B_field[:,0]**2 + B_field[:,1]**2 + B_field[:,2]**2)**0.5
    B_field_mag = pd.Series(dict(zip(time, B_field_mag)))

    fig = plt.figure()
    ax = fig.gca()
    plt.style.use('seaborn-v0_8-whitegrid')

    B_field_mag.plot(ax=ax, label="|B|", color="midnightblue")
    B_field_x.plot(ax=ax, label="Bx", color="royalblue")
    B_field_y.plot(ax=ax, label="By", color="orange")
    B_field_z.plot(ax=ax, label="Bz", color="deeppink")

    plt.legend()

    # plot radial distance of juice wrt callisto
    radial.plot(ax=ax, label="distance", color="k", secondary_y=True)
    ax.right_ax.set_ylabel("Radial distance [R_C]")
    plt.legend(loc="lower right")

    # add a vertical line at CA
    ax.axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)

    # Format the time on the x-axis to include minutes
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M'))

    ax.set_ylabel("B-field [nT]")
    ax.set_title(field_type + " field during Flyby %s" % (flyby_n))
    plt.show()



def plot_time_evolution_Gal(B_field, orbit_cphio, orbit_CA, flyby_n, field_type):

    # convert time into Timestamp format
    J2000 = datetime(2000,1,1,12) # difference between J2000 and UTC
    time = orbit_cphio[0]

    time = [Timestamp((J2000 + timedelta(seconds=timestamp)).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time]
    time_CA = Timestamp((J2000 + timedelta(seconds=orbit_CA[0])).strftime('%Y-%m-%d %H:%M:%S'))

    # radial distance of juice wrt callisto
    radial = pd.Series(dict(zip(time, orbit_cphio[4] / R_C)))

    B_field_x = pd.Series(dict(zip(time, B_field[:,0])))
    B_field_y = pd.Series(dict(zip(time, B_field[:,1])))
    B_field_z = pd.Series(dict(zip(time, B_field[:,2])))

    B_field_mag = (B_field[:,0]**2 + B_field[:,1]**2 + B_field[:,2]**2)**0.5
    B_field_mag = pd.Series(dict(zip(time, B_field_mag)))

    fig = plt.figure()
    ax = fig.gca()
    plt.style.use('seaborn-v0_8-whitegrid')

    B_field_mag.plot(ax=ax, label="|B|", color="midnightblue")
    B_field_x.plot(ax=ax, label="Bx", color="royalblue")
    B_field_y.plot(ax=ax, label="By", color="orange")
    B_field_z.plot(ax=ax, label="Bz", color="deeppink")

    plt.legend()

    # plot radial distance of juice wrt callisto
    radial.plot(ax=ax, label="distance", color="k", secondary_y=True)
    ax.right_ax.set_ylabel("Radial distance [R_C]")
    plt.legend(loc="lower right")

    # add a vertical line at CA
    ax.axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)

    ax.set_ylabel("B-field [nT]")
    ax.set_title(field_type + " field during Flyby %s" % (flyby_n))
    plt.show()