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
p_directory = directory + "/previous/"
n_directory = directory+ "/next/"

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

rinex_p = glob.glob(p_directory+"*.Cmn")
std_p = glob.glob(p_directory+"*.Std")
rinex_n = glob.glob(n_directory+"*.Cmn")
std_n = glob.glob(n_directory+"*.Std")

load_dir_p = [open(rinex_p[i], "r") for i in range(len(rinex_p))]
load_std_p = [Table.read(std_p[i], format="ascii") for i in range(len(std_p))]
load_dir_n = [open(rinex_n[i], "r") for i in range(len(rinex_n))]
load_std_n = [Table.read(std_n[i], format="ascii") for i in range(len(std_n))]

# Plot vTEC map
fig = plt.figure()
ax = fig.add_subplot(3, 2, 3, adjustable="box", aspect="equal")
ax1 = fig.add_subplot(3, 2, 4, adjustable="box")
axp = fig.add_subplot(3, 2, 1, adjustable="box", aspect="equal")
axp1 = fig.add_subplot(3, 2, 2, adjustable="box")
axn = fig.add_subplot(3, 2, 5, adjustable="box", aspect="equal")
axn1 = fig.add_subplot(3, 2, 6, adjustable="box")


# Load and plot event position and start time

load_meteor_pos = Table.read("meteors_database.tab", format="ascii")
meteor_mask = load_meteor_pos["Fecha"] == date
ax.plot(load_meteor_pos["Longitud"][meteor_mask], load_meteor_pos["Latitud"][meteor_mask], "mo")
#ax.annotate("Event", (load_meteor_pos["Longitud"][meteor_mask], load_meteor_pos["Latitud"][meteor_mask]),
#			textcoords="offset points", color="w", xytext=(10, 10), ha="center", bbox=dict(boxstyle="round", pad=0.5, fc="r", alpha=0.7))

t0_meteor_1 = load_meteor_pos["T_0 (GLM-16)"][meteor_mask]
t0_meteor_2 = load_meteor_pos["T_0 (GLM-17)"][meteor_mask]

if t0_meteor_1 == "N/A":
    t0_m1_h, t0_m1_m, t0_m1_s = (np.nan, np.nan, np.nan)
    t0_m2_h, t0_m2_m, t0_m2_s = t0_meteor_2[0].split(":")
elif t0_meteor_2 == "N/A":
    t0_m2_h, t0_m2_m, t0_m2_s = (np.nan, np.nan, np.nan)
    t0_m1_h, t0_m1_m, t0_m1_s = t0_meteor_1[0].split(":")
else: # convert start time from string to float (in hours)
    t0_m1_h, t0_m1_m, t0_m1_s = t0_meteor_1[0].split(":")
    t0_m2_h, t0_m2_m, t0_m2_s = t0_meteor_2[0].split(":")

t0_m1 = float(t0_m1_h) + float(t0_m1_m)/60. + float(t0_m1_s)/3600.
t0_m2 = float(t0_m2_h) + float(t0_m2_m)/60. + float(t0_m2_s)/3600.

# Load and plot RINEX data

for f, g, fp, gp, fn, gn in zip(load_dirs, load_std, load_dir_p, load_std_p, load_dir_n, load_std_n):
    header = f.readline()
    header_p = fp.readline()
    header_n = fn.readline()
    h1, h2 = header.split(",")
    station = h2.split("\\")[-1][0:4]
    blank = f.readline()
    blank = fp.readline()
    blank = fn.readline()
    s_coords = f.readline()
    s_coords_p = fp.readline()
    s_coords_n = fn.readline()
    s_latitude, s_longitude, s_altitude = s_coords.split()
    blank = f.readline()
    blank = fp.readline()
    blank = fn.readline()
    data  = f.readlines()
    data_p = fp.readlines()
    data_n = fn.readlines()
    obs_tab = Table.read(data, format="ascii")
    obs_tab_p = Table.read(data_p, format="ascii")
    obs_tab_n = Table.read(data_n, format="ascii")
    std_time = g["col1"]
    std_time_p = gp["col1"]
    std_time_n = gn["col1"]
    std_TEC = g["col2"]
    std_TEC_p = gp["col2"]
    std_TEC_n = gn["col2"]
    for i in range(len(obs_tab["Vtec"])): # Replace "-" into NaN since there is no data
        if obs_tab["Vtec"][i] == "-":
            obs_tab["Vtec"][i] = np.nan
    for i in range(len(obs_tab_p["Vtec"])):
        if obs_tab_p["Vtec"][i] == "-":
            obs_tab_p["Vtec"][i] = np.nan
    for i in range(len(obs_tab_n["Vtec"])):
        if obs_tab_n["Vtec"][i] == "-":
            obs_tab_n["Vtec"][i] = np.nan

    for i in range(len(std_TEC)):
        if std_TEC[i] == "-":
            std_TEC[i]=np.nan
    for i in range(len(std_TEC_p)):
        if std_TEC_p[i] == "-":
            std_TEC_p[i]=np.nan
    for i in range(len(std_TEC_n)):
        if std_TEC_n[i] == "-":
            std_TEC_n[i]=np.nan

    mean_TEC_int = interp1d(std_time, std_TEC)
    mean_TEC_int_p = interp1d(std_time_p, std_TEC_p)
    mean_TEC_int_n = interp1d(std_time_n, std_TEC_p)
    cmn_time = obs_tab["Time"]
    cmn_time_p = obs_tab_p["Time"]
    cmn_time_n = obs_tab_n["Time"]
    mask = cmn_time < 0
    mask_p = cmn_time_p < 0
    mask_n = cmn_time_n < 0
    cmn_time[mask] = cmn_time[mask] + 24.
    cmn_time_p[mask_p] = cmn_time_p[mask_p] + 24.0
    cmn_time_n[mask_n] = cmn_time_n[mask_n] + 24.0
    mask2 = cmn_time < max(std_time)
    mask2_p = cmn_time_p < max(std_time_p)
    mask2_n = cmn_time_n < max(std_time_n)
    dTEC = obs_tab["Vtec"][mask2] - mean_TEC_int(cmn_time[mask2])
    dTEC_p = obs_tab_p["Vtec"][mask2_p] - mean_TEC_int_p(cmn_time_p[mask2_p])
    dTEC_n = obs_tab_n["Vtec"][mask2_n] - mean_TEC_int_n(cmn_time_n[mask2_n])
    norm = MidpointNormalize(midpoint=0)
#    ax.plot(float(s_longitude)-360, float(s_latitude), "r*")
#    ax.text(float(s_longitude)-360+3, float(s_latitude), station.upper(), c="w",
#			bbox=dict(boxstyle='round', pad=0.5, fc='blue', alpha=0.3))
    im=ax.scatter(obs_tab["Lon"][mask2]-360, obs_tab["Lat"][mask2], s=1, c=dTEC, cmap="viridis",alpha=0.6, norm=norm)
    im1=ax1.scatter(cmn_time[mask2], obs_tab["Lat"][mask2], s=1, c=dTEC, cmap="viridis", alpha=0.6, norm=norm)
    im_p = axp.scatter(obs_tab_p["Lon"][mask2_p]-360, obs_tab_p["Lat"][mask2_p], s=1, c=dTEC_p, cmap="viridis", alpha=0.6, norm=norm)
    im1_p = axp1.scatter(cmn_time_p[mask2_p], obs_tab_p["Lat"][mask2_p], s=1, c=dTEC_p, cmap="viridis", alpha=0.6, norm=norm)
    im_n = axn.scatter(obs_tab_n["Lon"][mask2_n]-360, obs_tab_n["Lat"][mask2_n], s=1, c=dTEC_n, cmap="viridis", alpha=0.6, norm=norm)
    im1_n = axn1.scatter(cmn_time_n[mask2_n], obs_tab_n["Lat"][mask2_n], s=1, c=dTEC_n, cmap="viridis", alpha=0.6, norm=norm)


# Plot bolide trajectory

GLM16_file = open(directory+"/GLM/GLM-16-data.csv")
GLM17_file = open(directory+"/GLM/GLM-17-data.csv")

for i in range(10): # skip unneeded data
    GLM16_file.readline()
    GLM17_file.readline()

GLM16_data = GLM16_file.readlines()
GLM17_data = GLM17_file.readlines()
GLM16_table = Table.read(GLM16_data, format="ascii")
GLM17_table = Table.read(GLM17_data, format="ascii")

f1_longitude, f1_latitude = GLM16_table["longitude"], GLM16_table["latitude"]
f2_longitude, f2_latitude = GLM17_table["longitude"], GLM17_table["latitude"]


fit_coord1 = np.polyfit(f1_longitude, f1_latitude, 1)#nomial.Polynomial.fit(f1_longitude, f1_latitude, 1)
fit_coord2 = np.polyfit(f2_longitude, f2_latitude, 1)# polynomial.Polynomial.fit(f2_longitude, f2_latitude, 1)

xfit1 = np.linspace(51*f1_longitude[0]-50*f1_longitude[-1], f1_longitude[-1])
xfit2 = np.linspace(51*f2_longitude[0]-50*f2_longitude[-1], f2_longitude[-1])
yfit1 = f1_latitude[-1] + fit_coord1[0]*(xfit1-f1_longitude[-1])
yfit2 = f2_latitude[-1] + fit_coord2[0]*(xfit2-f2_longitude[-1])
#xfit1, yfit1 = fit_coord1.linspace()
#xfit2, yfit2 = fit_coord2.linspace()
x_trajectory, y_trajectory = 0.5*(xfit1+xfit2), 0.5*(yfit1+yfit2)
ax.plot(x_trajectory, y_trajectory, "k", lw=2)
ax.plot(xfit1, yfit1, "r--")
ax.plot(xfit2, yfit2, "r--")

# Show the interval of time the event started

ax1.axvline(x=t0_m1, ls="--", c="k")
ax1.axvline(x=t0_m2, ls="--", c="k")
ax1.axvspan(min((t0_m1, t0_m2)), max((t0_m1, t0_m2)), alpha=0.5, color="red")


# Plot settings

#ax = plt.gca()
#ax.set_aspect('equal', adjustable='box')
#plt.legend()
cbar = fig.colorbar(im, ax=ax)
cbar_p = fig.colorbar(im_p, ax=axp)
cbar_n = fig.colorbar(im_n, ax=axn)
cbar.set_label("Delta vTEC (TECU)")
cbar_p.set_label("Delta vTEC (TECU)")
cbar_n.set_label("Delta vTEC (TECU)")
cbar1 = fig.colorbar(im1, ax=ax1)
cbar1_p = fig.colorbar(im1_p, ax=axp1)
cbar1_n = fig.colorbar(im1_n, ax=axn1)
cbar1.set_label("Delta vTEC (TECU)")
cbar1_p.set_label("Delta vTEC (TECU)")
cbar1_n.set_label("Delta vTEC (TECU)")
out_dir = "./vTEC-maps/"
axn.set_xlabel("Longitude (deg)")
ax.set_ylabel("Latitude (deg)")
axp.set_ylabel("Latitude (deg)")
axn.set_ylabel("Latitude (deg)")
axn1.set_xlabel("Time (UT)")
plt.suptitle(date+" vTEC map")
ax1.set_ylabel("Latitude (deg)")
axp1.set_ylabel("Latitude (deg)")
axn1.set_ylabel("Latitude (deg)")
axp.title.set_text("A. Previous day")
ax.title.set_text("B. Event date")
axn.title.set_text("C. Next day")
fig.tight_layout()
fig.set_size_inches(18, 12)
plt.savefig(out_dir+date+"-vTEC_map.pdf")
