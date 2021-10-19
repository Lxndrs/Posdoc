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
from matplotlib.ticker import FormatStrFormatter


# This is a variant of the plot_vTEC.py script which filters the vTEC data several hours before the event
# and serveral hours after the event (terminal input) and plots the PRNs which satisfy those conditions.
# We expect to check if any of the chosen stations detected something related to the meteor pass

parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
				 help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--formato", type=str, default="pdf", choices=("pdf", "png", "jpg"), 
				  help="Choose output format")

parser.add_argument("--log", action="store_true", help="Use logarithmic scale for vTEC")

parser.add_argument("--filtertime", type=float, default=1.0, help="Time filtered after the event")

parser.add_argument("--filterbegins", type=float, default=0.2, help="Time filtered before the event")
parser.add_argument("--data", type=str, default="vTEC", choices=("vTEC", "sTEC", "ROTI"), help="Data type to plot")

cmd_args = parser.parse_args()
date = cmd_args.date
formato = cmd_args.formato
log = cmd_args.log
filtertime = cmd_args.filtertime
filterbegins=cmd_args.filterbegins
datatype = cmd_args.data


directory = "./data/"+date

# set figure style
sns.set_style("whitegrid") 


# Load RINEX capabilities

rinex_files = glob.glob(directory+"/*.Cmn")
std_files = glob.glob(directory+"/*.Std")
load_dirs = [open(rinex_files[i], "r") for i in range(len(rinex_files))]
load_std = [Table.read(std_files[i], format="ascii") for i in range(len(std_files))]
#load_ROTI = glob.glob(directory+"/ROTI*")
ROTI_files=[Table.read(load_ROTI[i], format="ascii") for i in range(len(load_ROTI))]

load_back = glob.glob(directory+"/*TEC.tab")
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
ax = fig.add_subplot(1, 2, 1, adjustable="box", aspect="equal")
ax1 = fig.add_subplot(1, 2, 2, adjustable="box")


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

# Load and plot RINEX data

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

    for i in range(len(std_TEC)):
        if std_TEC[i] == "-":
            std_TEC[i]=np.nan

    mean_TEC_int = interp1d(std_time, std_TEC)
    cmn_time = obs_tab["Time"]
    mask = cmn_time < 0
    cmn_time[mask] = cmn_time[mask] + 24.
    mask2 = (cmn_time > t0_m1 -filterbegins) & (cmn_time < t0_m1 + filtertime) # Now the mask filters the data in the desired interval time
    subs_tab = Table.read(directory+"/"+stations_dict[station], format="ascii")
    subs_TEC = subs_tab["Mean vTEC"]
    latitude, longitude = obs_tab["Lat"], obs_tab["Lon"]
    for i in range(len(subs_TEC)):
        if subs_TEC[i] == "-":
            subs_TEC[i] = np.nan
    subs_TEC_int = interp1d(std_time, subs_TEC) # Maybe substract the mean TEC AND the background TEC is redundant
    if datatype == "vTEC":
        dTEC =  obs_tab["Vtec"][mask2] - subs_TEC_int(cmn_time[mask2])
    elif datatype == "sTEC":
        dTEC = obs_tab["Stec"][mask2] - subs_TEC_int(cmn_time[mask2])
    norm = MidpointNormalize(midpoint=0)
    if log ==True: 
        im1 = ax1.scatter(cmn_time[mask2], latitude[mask2], s=1, c=dTEC[mask2], cmap="viridis", alpha=0.6, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, base=10))
        im=ax.scatter(longitude[mask2]-360, longitude[mask2], s=1, c=dTEC[mask2], cmap="viridis",alpha=0.6, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, base=10))
    else:
        im1=ax1.scatter(cmn_time[mask2],latitude[mask2], s=1, c=dTEC[mask2], cmap="viridis", alpha=0.6, norm=norm)
        im=ax.scatter(longitude[mask2]-360, latitude[mask2], s=1, c=dTEC[mask2], cmap="viridis",alpha=0.6, norm=norm)

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

#newlabel = [0, 5, 10, 15, 20, 25] Deprecated

time_zone_dict = {"2019-05-23":-5, "2019-07-18":-5, "2019-08-10":-5, "2019-10-03":-5, "2019-10-09":-7, "2019-11-16":-6, "2019-11-17":-8, "2019-11-19":-7, "2019-11-26":-5, "2019-12-04":-7, "2019-12-15":-7, "2019-12-29":-8, "2020-01-03":-8, "2020-01-06":-7, "2020-01-15":-6, "2020-02-12":-6, "2020-03-03":-6, "2020-03-31":-7, "2020-04-08":-5, "2020-04-18":-6, "2020-04-20":-5, "2020-04-25":-7, "2020-04-28":-6, "2020-05-08":-5, "2020-07-15":-6, "2020-08-07":-6, "2020-09-13":-7, "2020-09-30":-7, "2020-11-16":-6, "2020-11-17":-6, "2020-12-19":-6, "2020-12-23":-7, "2020-12-29":-6, "2021-03-31":-6}



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



ax1.axvspan(max(0, sunset_p_UT-24.0), sunrise_UT, alpha=0.1, color="cyan")
if sunset_p_UT-24 > 0:
    ax1.axvspan(0, sunset_p_UT-24, alpha=0.1, color="yellow")
ax1.axvspan(sunrise_UT, min(sunset_UT, 24.0), alpha=0.1, color="yellow")
if sunset_UT < 24.0:
    ax1.axvspan(sunset_UT, 24.0, alpha=0.1, color="cyan")

label = "Detrended TEC for {} stations".format(len(rinex_files))
if log ==True:
    cbar = fig.colorbar(im, ax=ax, ticks=[-1e1, -1, -1e-1, 1e-1, 1, 1e1])
    cbar.ax.minorticks_on()
    cbar1 = fig.colorbar(im1, ax=ax1, ticks=[-1e1, -1, -1e-1, 1e-1, 1, 1e1])
    cbar1.ax.minorticks_on()
else:
    cbar = fig.colorbar(im, ax=ax)
    cbar1 = fig.colorbar(im1, ax=ax1)

cbar.set_label("Delta {} (TECU)".format(datatype))
cbar1.set_label("Delta {} (TECU)".format(datatype))
if datatype=="ROTI":
    cbar.set_label(r"ROTI ($TECU~s^{-1}$)")
    cbar1.set_label(r"ROTI ($TECU~s^{-1}$)")
out_dir = "./vTEC-maps/"

ax.set_ylabel("Latitude (deg)")
ax1.set_xlabel("Universal Time (hours)")
plt.suptitle(date+" {} map".format(datatype))
ax1.set_ylabel("Latitude (deg)")
ax.title.set_text(label + ". Event date")
ax1.set_xlim(t0_m1 -filterbegins, t0_m1 + filtertime) # Change creation of second axis after setting xlim of original axis
oldlabel = ax1.get_xticks()
local_time = np.array(oldlabel) + time_zone_dict[date]
for i in range(len(local_time)):
    if local_time[i] < 0.:
        local_time[i] = local_time[i] + 24.0
ax2 = ax1.twiny()
ax2.set_xticks(oldlabel)
ax2.set_xticklabels(local_time)
ax2.set_xlabel("Local Time (hours)")
ax2.set_xlim(ax1.get_xlim())
ax1.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
ax2.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
fig.set_size_inches(11, 9)
fig.tight_layout()

if log==True:
    plt.savefig(out_dir+date+"-{}_logmap_filtered+{}-{}.".format(datatype, filtertime, filterbegins)+formato)
else:
    plt.savefig(out_dir+date+"-{}_map_filtered+{}-{}.".format(datatype, filtertime, filterbegins)+formato)
