import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp

from load_data import *
from induced_field import *
from jupiter_field import *


juice_wrt_callisto_cphio = get_spice_data("juice", "callisto", "cphio", "J")
callisto_wrt_jupiter_JSO = get_spice_data("juice", "jupiter", "jupsunorb", "J")

orbit4_cphio = juice_wrt_callisto_cphio["orbit4"]
orbit4_JSO = callisto_wrt_jupiter_JSO["orbit4"]

time = orbit4_cphio[0]
time = [Timestamp(datetime.utcfromtimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time]

B_external = Bext_full(orbit4_JSO)
B_induced = B_induced_infinite(orbit4_cphio, B_external, R_C, 100)

B_ind_x = pd.Series(dict(zip(time, B_induced[:,0])))
B_ind_y = pd.Series(dict(zip(time, B_induced[:,1])))
B_ind_z = pd.Series(dict(zip(time, B_induced[:,2])))

B_ext_x = pd.Series(dict(zip(time, B_external[:,0])))
B_ext_y = pd.Series(dict(zip(time, B_external[:,1])))
B_ext_z = pd.Series(dict(zip(time, B_external[:,2])))


fig = plt.figure()
ax = fig.gca()
plt.style.use('seaborn-v0_8-whitegrid')

B_ind_x.plot(ax=ax)
B_ind_y.plot(ax=ax)
B_ind_z.plot(ax=ax)

B_ext_x.plot(ax=ax)
B_ext_y.plot(ax=ax)
B_ext_z.plot(ax=ax)

plt.show()

