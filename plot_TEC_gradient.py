# Mexico map plotter
# The main idea of this program was taken from 
# https://towardsdatascience.com/mapping-with-matplotlib-pandas-geopandas-and-basemap-in-python-d11b57ab5dac
# By Ashwani Dhankhar 
# And the shape file for Mexico from CONABIO
# http://www.conabio.gob.mx/informacion/metadata/gis/destdv250k_2gw.xml?_xsl=/db/meadata/xsl/fgdc_html.xsl&_indent=no

# This is a clone of the vTEC_maps.py script except in this case we will compute the gradient of TEC instead of substracting the actual TEC from 
# some backgrond (or not?)

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

parser.add_argument("--formato", type=str, default="pdf", choices=("pdf", "png", "jpg"), 
                                help="Choose output format")

parser.add_argument("--log", action="store_true", help="Use logarithmic scale for vTEC")

parser.add_argument("--substract", action="store_true", help="substract the median of 27 previous days")
parser.add_argument("--data", type=str, default="vTEC", choices=("vTEC", "sTEC"))

cmd_args = parser.parse_args()
date = cmd_args.date
formato = cmd_args.formato
log = cmd_args.log
substract = cmd_args.substract
datatype = cmd_args.data
directory = "./data/"+date
p_directory = directory + "/previous/"
n_directory = directory+ "/next/"

# set figure style
sns.set_style("whitegrid") 


# Load RINEX capabilities

rinex_files = glob.glob(directory+"/*.Cmn")
#std_files = glob.glob(directory+"/*.Std")
load_dirs = [open(rinex_files[i], "r") for i in range(len(rinex_files))]
#load_std = [Table.read(std_files[i], format="ascii") for i in range(len(std_files))]


rinex_p = glob.glob(p_directory+"*.Cmn")
#std_p = glob.glob(p_directory+"*.Std")
rinex_n = glob.glob(n_directory+"*.Cmn")
#std_n = glob.glob(n_directory+"*.Std")

load_dir_p = [open(rinex_p[i], "r") for i in range(len(rinex_p))]
#load_std_p = [Table.read(std_p[i], format="ascii") for i in range(len(std_p))]
load_dir_n = [open(rinex_n[i], "r") for i in range(len(rinex_n))]
#load_std_n = [Table.read(std_n[i], format="ascii") for i in range(len(std_n))]

if substract ==True:
    load_back = glob.glob(directory+"/*.tab")
    stations_names = []
    stations_files = []
    for l in load_back:
        home, dfolder, fecha, tabfile = l.split("/")
        s_name = tabfile.split("-")[0]
        stations_names.append(s_name)
        stations_files.append(tabfile)
    stations_dict = dict(zip(stations_names, stations_files))

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

t0_meteor_1 = load_meteor_pos["T_0 (GLM-16)"][meteor_mask]
t0_meteor_2 = load_meteor_pos["T_0 (GLM-17)"][meteor_mask]

t0_m1_h, t0_m1_m, t0_m1_s = t0_meteor_1[0].split(":")
t0_m2_h, t0_m2_m, t0_m2_s = t0_meteor_2[0].split(":")

t0_m1 = float(t0_m1_h) + float(t0_m1_m)/60. + float(t0_m1_s)/3600.
t0_m2 = float(t0_m2_h) + float(t0_m2_m)/60. + float(t0_m2_s)/3600.
#stations = []
# Load and plot RINEX data

for f, fp, fn in zip(load_dirs, load_dir_p, load_dir_n):
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
#    std_time = g["col1"]
#    std_time_p = gp["col1"]
#    std_time_n = gn["col1"]
#    std_TEC = g["col2"]
#    std_TEC_p = gp["col2"]
#    std_TEC_n = gn["col2"]
    for i in range(len(obs_tab["Vtec"])): # Replace "-" into NaN since there is no data
        if obs_tab["Vtec"][i] == "-":
            obs_tab["Vtec"][i] = np.nan
    for i in range(len(obs_tab_p["Vtec"])):
        if obs_tab_p["Vtec"][i] == "-":
            obs_tab_p["Vtec"][i] = np.nan
    for i in range(len(obs_tab_n["Vtec"])):
        if obs_tab_n["Vtec"][i] == "-":
            obs_tab_n["Vtec"][i] = np.nan

#    for i in range(len(std_TEC)):
#        if std_TEC[i] == "-":
#            std_TEC[i]=np.nan
#    for i in range(len(std_TEC_p)):
#        if std_TEC_p[i] == "-":
#            std_TEC_p[i]=np.nan
#    for i in range(len(std_TEC_n)):
#        if std_TEC_n[i] == "-":
#            std_TEC_n[i]=np.nan

#    mean_TEC_int = interp1d(std_time, std_TEC)
#    mean_TEC_int_p = interp1d(std_time_p, std_TEC_p)
#    mean_TEC_int_n = interp1d(std_time_n, std_TEC_p)
    cmn_time = obs_tab["Time"]
    cmn_time_p = obs_tab_p["Time"]
    cmn_time_n = obs_tab_n["Time"]
    mask = cmn_time < 0
    mask_p = cmn_time_p < 0
    mask_n = cmn_time_n < 0
    cmn_time[mask] = cmn_time[mask] + 24.
    cmn_time_p[mask_p] = cmn_time_p[mask_p] + 24.0
    cmn_time_n[mask_n] = cmn_time_n[mask_n] + 24.0
#    mask2 = cmn_time < max(std_time)
#    mask2_p = cmn_time_p < max(std_time_p)
#    mask2_n = cmn_time_n < max(std_time_n)
    if substract==True:
        subs_tab = Table.read(directory+"/"+stations_dict[station], format="ascii")
        subs_TEC = subs_tab["Mean vTEC"]
        std_time = subs_tab["Time (UT)"]
        mask2 = cmn_time < max(std_time)
        mask2_p = cmn_time_p < max(std_time)
        mask2_n = cmn_time_n < max(std_time)
        for i in range(len(subs_TEC)):
            if subs_TEC[i] == "-":
                subs_TEC[i] = np.nan
        subs_TEC_int = interp1d(std_time, subs_TEC) # Maybe substract the mean TEC AND the background TEC is redundant
        time=cmn_time[mask2]
        time_p = cmn_time_p[mask2_p]
        time_n = cmn_time_n[mask2_n]
        latitud = obs_tab["Lat"][mask2]
        latitud_p = obs_tab_p["Lat"][mask2_p]
        latitud_n =obs_tab_n["Lat"][mask2_n]
        longitud = obs_tab["Lon"][mask2]      
        longitud_p = obs_tab_p["Lon"][mask2_p]
        longitud_n =obs_tab_n["Lon"][mask2_n] 
        if datatype == "vTEC":
            dTEC =  np.gradient(obs_tab["Vtec"][mask2] - subs_TEC_int(cmn_time[mask2]))/np.gradient(cmn_time[mask2])
            dTEC_p = np.gradient(obs_tab_p["Vtec"][mask2_p]-subs_TEC_int(cmn_time_p[mask2_p]))/np.gradient(cmn_time_p[mask2_p])
            dTEC_n = np.gradient(obs_tab_n["Vtec"][mask2_n]-subs_TEC_int(cmn_time_n[mask2_n]))/np.gradient(cmn_time_n[mask2_n])
        elif datatype == "sTEC":
            dTEC =  np.gradient(obs_tab["Stec"][mask2]-subs_TEC_int(cmn_time[mask2]))/np.gradient(cmn_time[mask2]) 
            dTEC_p = np.gradient(obs_tab_p["Stec"][mask2_p]-subs_TEC_int(cmn_time[mask2_p]))/np.gradient(cmn_time_p[mask2_p]) 
            dTEC_n = np.gradient(obs_tab_n["Stec"][mask2_n]-subs_TEC_int(cmn_time[mask2_n]))/np.gradient(cmn_time_n[mask2_n]) 
    else:
        time=cmn_time
        time_p = cmn_time_p
        time_n = cmn_time_n
        latitud = obs_tab["Lat"]
        latitud_p = obs_tab_p["Lat"]
        latitud_n =obs_tab_n["Lat"]
        longitud = obs_tab["Lon"]
        longitud_p = obs_tab_p["Lon"]
        longitud_n =obs_tab_n["Lon"]
        if datatype == "vTEC":
            dTEC = np.gradient(obs_tab["Vtec"])/np.gradient(cmn_time)
            dTEC_p = np.gradient(obs_tab_p["Vtec"])/np.gradient(cmn_time_p)
            dTEC_n = np.gradient(obs_tab_n["Vtec"])/np.gradient(cmn_time_n)
        elif datatype == "sTEC":
            dTEC = np.gradient(obs_tab["Vtec"])/np.gradient(cmn_time)
            dTEC_p = np.gradient(obs_tab_p["Vtec"])/np.gradient(cmn_time_p)
            dTEC_n = np.gradient(obs_tab_n["Vtec"])/np.gradient(cmn_time_n)
    if log ==False:
        norm = MidpointNormalize(midpoint=0)
    else:
        norm = colors.SymLogNorm(linthresh=0.03, linscale=0.03, base=10)
    im1 = ax1.scatter(time, latitud, s=1, c=dTEC, cmap="viridis", alpha=0.6, norm=norm)
    im1_p = axp1.scatter(time_p, latitud_p, s=1, c=dTEC_p, cmap="viridis", alpha=0.6, norm=norm)
    im1_n = axn1.scatter(time_n, latitud_n, s=1, c=dTEC_n, cmap="viridis", alpha=0.6, norm=norm)
    im = ax.scatter(longitud-360, latitud, s=1, c=dTEC, cmap="viridis",alpha=0.6, norm=norm)
    im_p = axp.scatter(longitud_p-360, latitud_p, s=1, c=dTEC_p, cmap="viridis", alpha=0.6, norm=norm)
    im_n = axn.scatter(longitud_n-360, latitud_n, s=1, c=dTEC_n, cmap="viridis", alpha=0.6, norm=norm)



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


fit_coord1 = np.polyfit(f1_longitude, f1_latitude, 1)
fit_coord2 = np.polyfit(f2_longitude, f2_latitude, 1)


poly1 = np.poly1d(fit_coord1)
poly2 = np.poly1d(fit_coord2)
if f1_longitude[-1] > f1_longitude[0]:
    step1 = -20
else:
    step1 = 20

if f2_longitude[-1] > f2_longitude[0]:
    step2 = -20
else:
    step2 = 20

ax.plot([f1_longitude[0]+step1, f1_longitude[-1]], [poly1(f1_longitude[0]+step1), poly1(f1_longitude[-1])], "r--")
ax.plot([f2_longitude[0]+step2, f2_longitude[-1]], [poly2(f2_longitude[0]+step2), poly2(f2_longitude[-1])], "r--")

# Show the interval of time the event started

ax1.axvline(x=t0_m1, ls="--", c="k")
ax1.axvline(x=t0_m2, ls="--", c="k")
ax1.axvspan(min((t0_m1, t0_m2)), max((t0_m1, t0_m2)), alpha=0.5, color="red")

# add the second x-axis 

newlabel = [0, 5, 10, 15, 20, 25]
time_zone_dict = {"2019-05-23":-5, "2019-07-18":-5, "2019-08-10":-5, "2019-10-03":-5, "2019-10-09":-7, "2019-11-16":-6, "2019-11-17":-8, "2019-11-19":-7, "2019-11-26":-5, "2019-12-04":-7, "2019-12-15":-7, "2019-12-29":-8, "2020-01-03":-8, "2020-01-06":-7, "2020-01-15":-6, "2020-02-12":-6, "2020-03-03":-6, "2020-03-31":-7, "2020-04-08":-5, "2020-04-18":-6, "2020-04-20":-5, "2020-04-25":-7, "2020-04-28":-6, "2020-05-08":-5, "2020-07-15":-6, "2020-08-07":-6, "2020-09-13":-7, "2020-09-30":-7, "2020-11-16":-6, "2020-11-17":-6, "2020-12-19":-6, "2020-12-23":-7, "2020-12-29":-6, "2021-03-31":-6}

local_time = np.array(newlabel) + time_zone_dict[date]
for i in range(len(local_time)):
    if local_time[i] < 0.:
        local_time[i] = local_time[i] + 24.0
ax2 = ax1.twiny()
ax2.set_xticks(newlabel)
ax2.set_xticklabels(local_time)
ax2.set_xlim(ax1.get_xlim())

axp2 = axp1.twiny()
axp2.set_xticks(newlabel)
axp2.set_xticklabels(local_time)
axp2.set_xlabel("Local Time (hours)")
axp2.set_xlim(axp1.get_xlim())

axn2 = axn1.twiny()
axn2.set_xticks(newlabel)
axn2.set_xticklabels(local_time)
axn2.set_xlim(axn1.get_xlim())

# Daytime shaded in light yellow and night time in blue

# local sunrise and sunset dictionaries 

sunrise_p = {"2019-05-23":6+57./60, "2019-07-18":7+14./60, "2019-08-10":7+25./60, "2019-10-03":7+26/60., "2019-10-09":7+16./60, "2019-11-16": 7+18./60, "2019-11-17": 6+19/60., "2019-11-19": 7+ 6./60, "2019-11-26": 6+33/60., "2019-12-04": 7+17./60, "2019-12-15": 7+0./60, "2019-12-29": 6+46/60., "2020-01-03": 6+47./60, "2020-01-06": 7+ 7./60, "2020-01-15": 7+1./60, "2020-02-12": 6+41./60, "2020-03-03": 7+16./60, "2020-03-31": 6+16/60., "2020-04-08": 7+ 4./60, "2020-04-18": 6+37./60, "2020-04-20": 7+3./60, "2020-04-25": 5+44./60, "2020-04-28": 6+50./60, "2020-05-08": 6+42./60, "2020-07-15":6.5, "2020-08-07":6.5, "2020-09-13": 7+14./60, "2020-09-30":7+13./60, "2020-11-16": 6+52./60, "2020-11-17": 7+5./60, "2020-12-19":7+20./60, "2020-12-23": 7+ 1./60, "2020-12-29": 7+19/60., "2021-03-31": 6+ 8./60}


sunset_p = {"2019-05-23":20+24./60, "2019-07-18":20+45./60, "2019-08-10":20+23./60, "2019-10-03":19+17/60., "2019-10-09":19, "2019-11-16":17+52./60, "2019-11-17":16+45/60., "2019-11-19":18+14./60, "2019-11-26":17+19/60., "2019-12-04":17.5, "2019-12-15":20+6./60, "2019-12-29":16+49/60., "2020-01-03":16+53./60, "2020-01-06":17+15./60, "2020-01-15":18+5./60, "2020-02-12":18+10./60, "2020-03-03":19+2./60, "2020-03-31":18+40/60., "2020-04-08":19+42./60, "2020-04-18":19.5, "2020-04-20":20+2./60, "2020-04-25":18+58./60, "2020-04-28":19+47./60, "2020-05-08":19+33./60, "2020-07-15":20, "2020-08-07":19+49./60, "2020-09-13":19+35./60, "2020-09-30":19+8./60, "2020-11-16":17+59./60, "2020-11-17":18+4./60, "2020-12-19":18+5./60, "2020-12-23":17+38./60, "2020-12-29":18+21/60., "2021-03-31":18+23./60}

sunrise = {"2019-05-23":6+56./60, "2019-07-18":7+14./60, "2019-08-10":7+26./60, "2019-10-03":7+27/60., "2019-10-09":7+17./60, "2019-11-16": 7+19./60, "2019-11-17": 6+20/60., "2019-11-19": 7+ 6./60, "2019-11-26": 6+34/60., "2019-12-04": 7+18./60, "2019-12-15": 7+1./60, "2019-12-29": 6+46/60., "2020-01-03": 6+47./60, "2020-01-06": 7+ 7./60, "2020-01-15": 7+1./60, "2020-02-12": 6+41./60, "2020-03-03": 7+15./60, "2020-03-31": 6+15/60., "2020-04-08": 7+ 3./60, "2020-04-18": 6+36./60, "2020-04-20": 7+2./60, "2020-04-25": 5+43./60, "2020-04-28": 6+49./60, "2020-05-08": 6+42./60, "2020-07-15":6.5, "2020-08-07":6.5, "2020-09-13": 7+14./60, "2020-09-30":7+13./60, "2020-11-16": 6+53./60, "2020-11-17": 7+3./60, "2020-12-19":7+19./60, "2020-12-23": 7+ 2./60, "2020-12-29": 7+20/60., "2021-03-31": 6+ 7./60}

sunset = {"2019-05-23":20+25./60, "2019-07-18":20+45./60, "2019-08-10":20+22./60, "2019-10-03":19+16/60., "2019-10-09":18+59./60, "2019-11-16":17+52./60, "2019-11-17":16+45/60., "2019-11-19":18+14./60, "2019-11-26":17+19/60., "2019-12-04":17.5, "2019-12-15":20+5./60, "2019-12-29":16+50/60., "2020-01-03":16+53./60, "2020-01-06":17+15./60, "2020-01-15":18+5./60, "2020-02-12":18+10./60, "2020-03-03":19+2./60, "2020-03-31":18+40/60., "2020-04-08":19+43./60, "2020-04-18":19.5, "2020-04-20":20+3./60, "2020-04-25":18+59./60, "2020-04-28":19+48./60, "2020-05-08":19+34./60, "2020-07-15":20, "2020-08-07":19+49./60, "2020-09-13":19+34./60, "2020-09-30":19+7./60, "2020-11-16":17+59./60, "2020-11-17":18+4./60, "2020-12-19":18+6./60, "2020-12-23":17+39./60, "2020-12-29":18+22/60., "2021-03-31":18+24./60}

sunrise_n= {"2019-05-23":6+56./60, "2019-07-18":7+15./60, "2019-08-10":7+26./60, "2019-10-03":7+27/60., "2019-10-09":7+17./60, "2019-11-16": 7+20./60, "2019-11-17": 6+21/60., "2019-11-19": 7+ 6./60, "2019-11-26": 6+35/60., "2019-12-04": 7+19./60, "2019-12-15": 7+1./60, "2019-12-29": 6+46/60., "2020-01-03": 6+48./60, "2020-01-06": 7+ 7./60, "2020-01-15": 7+8./60, "2020-02-12": 6+40./60, "2020-03-03": 7+14./60, "2020-03-31": 6+13/60., "2020-04-08": 7+ 1./60, "2020-04-18": 6+35./60, "2020-04-20": 7+1./60, "2020-04-25": 5+42./60, "2020-04-28": 6+49./60, "2020-05-08": 6+41./60, "2020-07-15":6.5, "2020-08-07":6.5, "2020-09-13": 7+15./60, "2020-09-30":7+14./60, "2020-11-16": 6+53./60, "2020-11-17": 7+6./60, "2020-12-19":7+20./60, "2020-12-23": 7+ 2./60, "2020-12-29": 7+20/60., "2021-03-31": 6+ 6./60}

sunset_n = {"2019-05-23":20+25./60, "2019-07-18":20+44./60, "2019-08-10":20+21./60, "2019-10-03":19+15/60., "2019-10-09":18+58./60, "2019-11-16":17+52./60, "2019-11-17":16+44/60., "2019-11-19":18+13./60, "2019-11-26":17+19/60., "2019-12-04":17.5, "2019-12-15":20+4./60, "2019-12-29":16+51/60., "2020-01-03":16+54./60, "2020-01-06":17+16./60, "2020-01-15":18+5./60, "2020-02-12":18+10./60, "2020-03-03":19+3./60, "2020-03-31":18+41/60., "2020-04-08":19+44./60, "2020-04-18":19.5, "2020-04-20":20+3./60, "2020-04-25":19, "2020-04-28":19+48./60, "2020-05-08":19+34./60, "2020-07-15":20, "2020-08-07":19+48./60, "2020-09-13":19+33./60, "2020-09-30":19+6./60, "2020-11-16":17+59./60, "2020-11-17":18+4./60, "2020-12-19":18+7./60, "2020-12-23":17+39./60, "2020-12-29":18+22/60., "2021-03-31":18+24./60}

sunrise_p_UT = sunrise_p[date] - time_zone_dict[date]
sunset_p_UT = sunset_p[date] - time_zone_dict[date]
sunrise_UT = sunrise[date] - time_zone_dict[date]
sunset_UT = sunset[date] - time_zone_dict[date]
sunrise_n_UT = sunrise_n[date] - time_zone_dict[date]
sunset_n_UT = sunset_n[date] - time_zone_dict[date]


axp1.axvspan(0, sunrise_p_UT, alpha=0.1, color="cyan")
axp1.axvspan(sunrise_p_UT, min(sunset_p_UT, 24.0), alpha=0.1, color="yellow")
if sunset_p_UT < 24.0:
    axp1.axvspan(sunset_p_UT, 24.0, alpha=0.1, color="cyan")
ax1.axvspan(max(0, sunset_p_UT-24.0), sunrise_UT, alpha=0.1, color="cyan")
if sunset_p_UT-24 > 0:
    ax1.axvspan(0, sunset_p_UT-24, alpha=0.1, color="yellow")
ax1.axvspan(sunrise_UT, min(sunset_UT, 24.0), alpha=0.1, color="yellow")
if sunset_UT < 24.0:
    ax1.axvspan(sunset_UT, 24.0, alpha=0.1, color="cyan")
axn1.axvspan(max(0, sunset_UT-24.0), sunrise_n_UT, alpha=0.1, color="cyan")
if sunset_UT - 24.0 >0:
    axn1.axvspan(0, sunset_UT-24, alpha=0.1, color="yellow")
axn1.axvspan(sunrise_n_UT, 24.0, alpha=0.1, color="yellow")

# Plot settings


label = "Rate of {} for {} stations".format(datatype, len(rinex_files))
if log ==True:
    cbar = fig.colorbar(im, ax=ax, ticks=[-1e1, -1, -1e-1, 1e-1, 1, 1e1])
    cbar.ax.minorticks_on()
    cbar_p = fig.colorbar(im_p, ax=axp, ticks=[-1e1, -1, -1e-1, 1e-1, 1, 1e1])
    cbar_p.ax.minorticks_on()
    cbar_n = fig.colorbar(im_n, ax=axn, ticks=[-1e1, -1, -1e-1, 1e-1, 1, 1e1])
    cbar_n.ax.minorticks_on()
    cbar1 = fig.colorbar(im1, ax=ax1, ticks=[-1e1, -1, -1e-1, 1e-1, 1, 1e1])
    cbar1.ax.minorticks_on()
    cbar1_p = fig.colorbar(im1_p, ax=axp1, ticks=[-1e1, -1, -1e-1, 1e-1, 1, 1e1])
    cbar1_p.ax.minorticks_on()
    cbar1_n = fig.colorbar(im1_n, ax=axn1, ticks=[-1e1, -1, -1e-1, 1e-1, 1, 1e1])
    cbar1_n.ax.minorticks_on()
else:
    cbar = fig.colorbar(im, ax=ax)
    cbar_p = fig.colorbar(im_p, ax=axp)
    cbar_n = fig.colorbar(im_n, ax=axn)
    cbar1 = fig.colorbar(im1, ax=ax1)
    cbar1_p = fig.colorbar(im1_p, ax=axp1)
    cbar1_n = fig.colorbar(im1_n, ax=axn1)

cbar.set_label("Rate of {} (TECU/s)".format(datatype))
cbar_p.set_label("Rate of {} (TECU/s)".format(datatype))
cbar_n.set_label("Rate of {} (TECU/s)".format(datatype))
cbar1.set_label("Rate of {} (TECU/s)".format(datatype))
cbar1_p.set_label("Rate of {} (TECU/s)".format(datatype))
cbar1_n.set_label("Rate of {} (TECU/s)".format(datatype))
out_dir = "./rTEC-maps/"
axn.set_xlabel("Longitude (deg)")
ax.set_ylabel("Latitude (deg)")
axp.set_ylabel("Latitude (deg)")
axn.set_ylabel("Latitude (deg)")
axn1.set_xlabel("Universal Time (hours)")
plt.suptitle(date+"Rate of {} map".format(datatype))
ax1.set_ylabel("Latitude (deg)")
axp1.set_ylabel("Latitude (deg)")
axn1.set_ylabel("Latitude (deg)")
axp.title.set_text(label + ". Previous day")
axp1.title.set_text(label)
ax.title.set_text(label + ". Event date")
#ax1.title.set_text(label)
axn.title.set_text(label+". Next day")
#axn1.title.set_text(label)
fig.set_size_inches(22, 18)
#fig.tight_layout()

if substract == True:
    if log ==True:
        plt.savefig(out_dir+date+"-Rate of {}_logmap_minus_background.".format(datatype)+formato)
    else:
        plt.savefig(out_dir+date+"-Rate of {}_map_minus_background.".format(datatype)+formato)
else:
    if log==True:
        plt.savefig(out_dir+date+"-Rate of {}_logmap.".format(datatype)+formato)
    else:
        plt.savefig(out_dir+date+"-Rate of {}_map.".format(datatype)+formato)
