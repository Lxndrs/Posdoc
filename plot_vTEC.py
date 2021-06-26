# Mexico map plotter
# The main idea of this program was taken from 
# https://towardsdatascience.com/mapping-with-matplotlib-pandas-geopandas-and-basemap-in-python-d11b57ab5dac
# By Ashwani Dhankhar 
# And the shape file for Mexico from CONABIO
# http://www.conabio.gob.mx/informacion/metadata/gis/destdv250k_2gw.xml?_xsl=/db/meadata/xsl/fgdc_html.xsl&_indent=no

import seaborn as sns
import numpy as np
import pandas as pd
import shapefile as shp
import matplotlib.pyplot as plt
from plotfullmap import plot_map
import argparse
from astropy.table import Table
import glob
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.interpolate import interp1d
from midpoint import MidpointNormalize



parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
			 help='Choose date. Format: yyyy-mm-dd')



cmd_args = parser.parse_args()
date = cmd_args.date


directory = "./data/"+date

# set figure style
sns.set_style("whitegrid") 
#sns.mpl.rc("figure", figsize=(10,6))

# Read shape file of Mexico map (deprecated)
#sf = shp.Reader("map.shp")
#plot_map(sf)


# Load RINEX capabilities

rinex_files = glob.glob(directory+"/*.Cmn")
std_files = glob.glob(directory+"/*.Std")
load_dirs = [open(rinex_files[i], "r") for i in range(len(rinex_files))]
load_std = [Table.read(std_files[i], format="ascii") for i in range(len(std_files))]

# Plot vTEC map
fig = plt.figure()
ax = fig.add_subplot(1, 2, 1, adjustable="box", aspect="equal")
ax1 = fig.add_subplot(1, 2, 2, adjustable="box")

# Load and plot event position

load_meteor_pos = Table.read("meteors_database.tab", format="ascii")
meteor_mask = load_meteor_pos["Fecha"] == date
ax.plot(load_meteor_pos["Longitud"][meteor_mask], load_meteor_pos["Latitud"][meteor_mask], "mo")
ax.annotate("Event", (load_meteor_pos["Longitud"][meteor_mask], load_meteor_pos["Latitud"][meteor_mask]),
		  textcoords="offset points", color="w", xytext=(10, 10), ha="center", bbox=dict(boxstyle="round", pad=0.5, fc="r", alpha=0.7))



for f, g in zip(load_dirs, load_std):
    header = f.readline()
    h1, h2 = header.split(",")
    station = h2.split("\\")[-1][0:4]
    blank = f.readline()
    s_coords = f.readline()
    s_latitude, s_longitude, s_altitude = s_coords.split()
    blank = f.readline()
    data  = f.readlines()
    obs_tab = Table.read(data, format="ascii")
    std_time = g["col1"]
    std_TEC = g["col2"]
    for i in range(len(obs_tab["Vtec"])): # Replace "-" into NaN since there is no data
        if obs_tab["Vtec"][i] == "-":
            obs_tab["Vtec"][i] = np.nan
    std_mask = std_TEC == "-" # Mask data with "-"
    mean_TEC_int = interp1d(std_time[~std_mask], std_TEC[~std_mask])
    cmn_time = obs_tab["Time"]
    mask = cmn_time < 0
    cmn_time[mask] = cmn_time[mask] + 24.
    mask2 = cmn_time < max(std_time)
    dTEC = obs_tab["Vtec"][mask2] - mean_TEC_int(cmn_time[mask2])
    norm = MidpointNormalize(midpoint=0)
    ax.plot(float(s_longitude)-360, float(s_latitude), "r*")
    ax.text(float(s_longitude)-360+3, float(s_latitude), station.upper(), c="w",
		  bbox=dict(boxstyle='round', pad=0.5, fc='blue', alpha=0.3))
    im=ax.scatter(obs_tab["Lon"][mask2]-360, obs_tab["Lat"][mask2], s=1, c=dTEC, cmap="viridis",alpha=0.8, norm=norm)
    im1=ax1.scatter(cmn_time[mask2], obs_tab["Lat"][mask2], s=1, c=dTEC, cmap="viridis", alpha=0.8, norm=norm)


# Plot settings

#ax = plt.gca()
#ax.set_aspect('equal', adjustable='box')
#plt.legend()
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Delta vTEC (TECU)")
cbar1 = fig.colorbar(im1, ax=ax1)
cbar1.set_label("Delta vTEC (TECU)")
out_dir = "./vTEC-maps/"
ax.set_xlabel("Longitude (deg)")
ax.set_ylabel("Latitude (deg)")
ax1.set_xlabel("Time (UT)")
plt.suptitle(date+" vTEC map")
ax1.set_ylabel("Laitude (deg)")
fig.tight_layout()
fig.set_size_inches(18, 9)
plt.savefig(out_dir+date+"-vTEC_map.pdf")
