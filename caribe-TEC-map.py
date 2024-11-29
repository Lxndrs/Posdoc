# Creation of the Caribbean map with the respective sTEC spatial distribution density ~~~~~~~~~~~~~

# Import all python relevant libraries for this work ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotfullmap import plot_map
from astropy.table import Table 
import glob
import matplotlib.cm as cm
from scipy.interpolate import griddata, interp1d
from midpoint import MidpointNormalize
import datetime
from scipy.integrate import quad
import argparse


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# In order to save time from compiling a lot of times the same information we will do some python 
# scripts and save outputsus to some files and read these files from here and use them repeteadly
# instead of compiling everything from zero each time.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract data from the command line

parser = argparse.ArgumentParser(description= """Command line arguments""")
parser.add_argument("--time", type=str, default="21:25:48", help="Choose time instant")
parser.add_argument("--TEC", type=str, default="vTEC", choices=("vTEC", "sTEC"), 
                    help="Choose type of TEC data")
cmd_args = parser.parse_args()
time_snap = cmd_args.time
TEC_choice = cmd_args.TEC
TEC_dict = {"vTEC":"TEC", "sTEC":"sTEC"}
hrs, mins, secs = time_snap.split(":")
time_snap_dec = float(hrs) + float(mins)/60. + float(secs)/3600.
sns.set_style("whitegrid")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load caribbean map 

car_map = Table.read("caribbean-map.tab", format="ascii")
plt.plot(car_map["Longitude"], car_map["Latitude"], "k.", ms=0.025)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load meteor trajectory 
# We want to plot the meteor trajectory **UNTIL** the current time. If this times occurs before
# passage, plots nothing, and if occurs after, then plot the whole trajectory
meteor_trajectory = Table.read("meteor-trajectory.tab", format="ascii")
trajectory_filter = meteor_trajectory["UT"] < time_snap_dec
xm = meteor_trajectory["Longitude (deg)"][trajectory_filter]
ym= meteor_trajectory["Latitude (deg)"][trajectory_filter]
plt.plot(xm, ym, "b--")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Load TEC table for time t
#
file_folder = "./data/2019-06-22/"
TEC_file = file_folder+"{}_snapshot_{}.tab".format(TEC_choice, time_snap)
TEC_table = Table.read(TEC_file, format="ascii")
x, y = TEC_table["Longitude"], TEC_table["Latitude"]
TEC = TEC_table[TEC_choice]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Plot the TEC data

## Set map boudaries and create grid
xmin, ymin = np.min(x), np.min(y)
xmax, ymax = np.max(x), np.max(y)
res = 0.1
nx = int((xmax-xmin)/res)*1j
ny = int((ymax-ymin)/res)*1j
grid_x, grid_y = np.mgrid[xmin:xmax:nx, ymin:ymax:ny]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
##Create coordinates array
#

coordinates = []
for X, Y in zip(x, y):
    coord_pair = [X, Y]
    coordinates.append(coord_pair)

grid0 = griddata(np.array(coordinates), TEC, (grid_x, grid_y), method="cubic")
norm = MidpointNormalize(midpoint=0)
plt.imshow(grid0.T, extent=(xmin, xmax, ymin, ymax), origin="lower", cmap="viridis", norm=norm)
cbar = plt.colorbar()
cbar.set_label(r"$\Delta$ {} (TECU)".format(TEC_choice))
plt.xlabel("Longitude (deg)")
plt.ylabel("Latitude (deg)")
plt.xlim(-80, -45)
plt.ylim(-5, 28)
plt.gcf().set_size_inches(10, 10)
outfile = TEC_file.replace(".tab", ".png")
plt.savefig(outfile)
