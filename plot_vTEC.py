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
from scipy.interpolate import interp1d



parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
		    help='Choose date. Format: yyyy-mm-dd')



cmd_args = parser.parse_args()
date = cmd_args.date


directory = "./data/"+date

# set figure style
sns.set_style("whitegrid") 
sns.mpl.rc("figure", figsize=(10,6))

# Read shape file of Mexico map
sf = shp.Reader("map.shp")
plot_map(sf)


# Load RINEX capabilities

rinex_files = glob.glob(directory+"/*.Cmn")
std_files = glob.glob(directory+"/*.Std")
load_dirs = [open(rinex_files[i], "r") for i in range(len(rinex_files))]
load_std = [Table.read(std_files[i], format="ascii") for i in range(len(std_files))]

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
    mean_TEC_int = interp1d(std_time, std_TEC)
    cmn_time = obs_tab["Time"]
    mask = cmn_time < 0
    cmn_time[mask] = cmn_time[mask] + 24.
    mask2 = cmn_time < max(std_time)
    dTEC = obs_tab["Vtec"][mask2] - mean_TEC_int(cmn_time[mask2])
    plt.plot(float(s_longitude)-360, float(s_latitude), "r*")
    plt.text(float(s_longitude)-360+0.5, float(s_latitude)-0.5, station.upper(),
	     bbox=dict(boxstyle='round', pad=0.5, fc='blue', alpha=1))

    plt.scatter(obs_tab["Lon"][mask2]-360, obs_tab["Lat"][mask2], s=1, c=dTEC, cmap="plasma",alpha=0.8)

# Plot settings

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.legend()
plt.colorbar()

plt.savefig(directory+"/"+date+"-GLM_map.pdf")
