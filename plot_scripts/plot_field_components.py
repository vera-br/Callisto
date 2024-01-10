# plot script for field components evolving in time

import numpy as np
import pandas as pd
from pandas import Timestamp
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator

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

    B_field_mag.plot(ax=ax, label="|B|", color="midnightblue")
    B_field_x.plot(ax=ax, label="Bx", color="royalblue")
    B_field_y.plot(ax=ax, label="By", color="orange")
    B_field_z.plot(ax=ax, label="Bz", color="deeppink")

    plt.legend(framealpha=1)

    # plot radial distance of juice wrt callisto
    radial.plot(ax=ax, label="distance", color="k", secondary_y=True)
    ax.right_ax.set_ylabel("Radial distance [R_C]")
    plt.legend(loc="lower right")

    # add a vertical line at CA
    ax.axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)

    ax.tick_params(axis='both', direction='in',top = True, right = True, which='major')
    ax.yaxis.set_minor_locator(MaxNLocator(50)) 
    #ax.xaxis.set_minor_locator(MaxNLocator(50))
    
    ax.grid(color='xkcd:dark blue',alpha =0.2)

    ax.set_ylabel("B-field [nT]")
    ax.set_title(field_type + " field during Flyby %s" % (flyby_n))
    plt.show()

def plot_compare_model_with_data(B_model, PDS_data, orbit_cphio, orbit_CA, title):
    """
    
    """

    # calculate B-field magnitudes
    B_mag = np.sqrt(PDS_data[1]**2 + PDS_data[2]**2 + PDS_data[3]**2)
    B_mag_model = np.sqrt(B_model[:,0]**2 + B_model[:,1]**2 + B_model[:,2]**2)

    # convert time from J2000 into Timestamp format
    J2000 = datetime(2000,1,1,12)

    time_CA = orbit_CA[0]
    time_CA = Timestamp((J2000 + timedelta(seconds=orbit_CA[0])).strftime('%Y-%m-%d %H:%M:%S'))

    time_model = orbit_cphio[0]
    time_model = [Timestamp((J2000 + timedelta(seconds=timestamp)).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time_model]

    time_Gal = PDS_data[0]
    time_Gal = [Timestamp((J2000 + timedelta(seconds=timestamp)).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time_Gal]

    # create series for plotting time evolution
    Bx_Gal = pd.Series(dict(zip(time_Gal, PDS_data[1])))
    By_Gal = pd.Series(dict(zip(time_Gal, PDS_data[2])))
    Bz_Gal = pd.Series(dict(zip(time_Gal, PDS_data[3])))
    Bmag_Gal = pd.Series(dict(zip(time_Gal, B_mag)))

    Bx_model = pd.Series(dict(zip(time_model, B_model[:, 0])))
    By_model = pd.Series(dict(zip(time_model, B_model[:, 1])))
    Bz_model = pd.Series(dict(zip(time_model, B_model[:, 2])))
    Bmag_model = pd.Series(dict(zip(time_model, B_mag_model)))

    # plot
    fig, ax = plt.subplots(2, 2, figsize=(12,5))

    Bx_Gal.plot(ax=ax[0,0], label="Data", color="skyblue")
    Bx_model.plot(ax=ax[0,0], label="Model", color="midnightblue")
    ax[0,0].set_title('Bx')

    By_Gal.plot(ax=ax[0,1], label="Data", color="skyblue")
    By_model.plot(ax=ax[0,1], label="Model", color="midnightblue")
    ax[0,1].set_title('By')

    Bz_Gal.plot(ax=ax[1,0], label="Data", color="skyblue")
    Bz_model.plot(ax=ax[1,0], label="Model", color="midnightblue")
    ax[1,0].set_title('Bz')

    Bmag_Gal.plot(ax=ax[1,1], label="Data", color="skyblue")
    Bmag_model.plot(ax=ax[1,1], label="Model", color="midnightblue")
    ax[1,1].set_title('|B|')

    ax[0,1].legend(loc="upper right", framealpha=1)

    # add a vertical line at CA
    ax[0,0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
    ax[0,1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
    ax[1,0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
    ax[1,1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)

    # set axes labels
    ax[0,0].set_ylabel("nT")
    ax[0,1].set_ylabel("nT")
    ax[1,0].set_ylabel("nT")
    ax[1,1].set_ylabel("nT")

    # set axes ticks
    ax[0,0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
    ax[0,0].yaxis.set_minor_locator(MaxNLocator(50))
    #ax[0,0].xaxis.set_minor_locator(MaxNLocator(50))

    ax[0,1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
    ax[0,1].yaxis.set_minor_locator(MaxNLocator(50))
    #ax[0,1].xaxis.set_minor_locator(MaxNLocator(50))

    ax[1,0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
    ax[1,0].yaxis.set_minor_locator(MaxNLocator(50)) 
    #ax[1,0].xaxis.set_minor_locator(MaxNLocator(50))

    ax[1,1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
    ax[1,1].yaxis.set_minor_locator(MaxNLocator(50))
    #ax[1,1].xaxis.set_minor_locator(MaxNLocator(50))

    # set grid
    ax[0,0].grid(color='xkcd:dark blue',alpha =0.2)
    ax[0,1].grid(color='xkcd:dark blue',alpha =0.2)
    ax[1,0].grid(color='xkcd:dark blue',alpha =0.2)
    ax[1,1].grid(color='xkcd:dark blue',alpha =0.2)

    fig.suptitle(title)
    plt.tight_layout()
    plt.show()