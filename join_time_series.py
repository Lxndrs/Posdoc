# This program has the purpose to attach two time series from adjacent dates from an
# specific station and an specific PRN and make a single time series.
# Starts as a clone of the time series grapher but I will start making changes

import pandas as pd
import argparse 
import numpy as np
import matplotlib.pyplot as plt
import glob
from statistics import mode
import seaborn as sns
parser = argparse.ArgumentParser(
	       description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
					  help='Choose date. Format: yyyy-mm-dd')


parser.add_argument("--previous", action="store_true", help="Join with previous day data")
parser.add_argument("--Next", action="store_true", help="Join with next day data")
parser.add_argument("--sets", type=str, help="Choose between sets of stations")
parser.add_argument("--sTEC", action="store_true", help="Use sTEC data")
parser.add_argument("--PRNs", type=str, default="1", help="Select PRNs which we need to join with another day")
parser.add_argument("--station", type=str, default="bara", help="Select station")

cmd_args = parser.parse_args()
date = cmd_args.date
previous = cmd_args.previous
Next = cmd_args.Next
stations_set = cmd_args.sets
sTEC = cmd_args.sTEC
set_folder = "/set"+stations_set
directory = "./data/"+date+set_folder
outdir = "./TEC_series/"+date+set_folder
PRNs = cmd_args.PRNs
station = cmd_args.station
if sTEC:
    infile = glob.glob(directory+"/*sTEC.csv")[0]
    title_label = "Detrended slant TEC series for individual satellites. Station {}".format(station.upper())
    g_extention = "/{}-joint-sTEC_series.pdf".format(station)
else:
    infile = glob.glob(directory+"/*-TEC.csv")[0]
    title_label = "Detrended vTEC series for individual satellites. Station {}".format(station.upper())
    g_extention = "/{}-joint-TEC_series.pdf".format(station)



# ****************************** Read data file ************************************************

csv_file = pd.read_csv(infile)
Date, time, residual = infile.split("_")
time_hours = float(time)
x_low = time_hours - 1.1
x_high = time_hours +200./60. +0.1

if previous:
    join_dir = directory+"/previous/"
    join_file = glob.glob(join_dir+"/-TEC.csv")[0]
    x_low = 22.0
    if sTEC:
        join_file = glob.glob(join_dir+"/*sTEC.csv")[0]
elif Next:
    join_dir = directory+"/next/"
    join_file = glob.glob(join_dir+"*-TEC.csv")[0]
    g_extention = "/{}-joint-TEC_series.pdf".format(station)
    x_high = 26.1
    if sTEC:
        join_file = glob.glob(join_dir+"*sTEC.csv")[0]
        g_extention = "/{}-joint-sTEC_series.pdf".format(station)
else:
    print("Error: must choose next or previous day")
    raise SystemExit
csv_join = pd.read_csv(join_file)
station_mask = csv_file["Station"] == station
PRN_set_str = PRNs.split(",")
PRN_set = [int(i) for i in PRN_set_str]
sns.set_style("whitegrid")
f = plt.figure()
n_subplots=len(PRN_set)
for i,p in enumerate(PRN_set):
    prn_mask = csv_file["PRN"] == p
    jprn_mask = csv_join["PRN"] ==p
    print("Station:{}, PRN:{}".format(station, p))
    time = csv_file["Time"][prn_mask & station_mask]
    jtime = csv_join["Time"][jprn_mask & station_mask]
    ylabel = "vTEC (TECU)"
    out_extention = "/{}-joint-TEC_series_PRN_{}.pdf".format(station, p)
    if sTEC:
        tec = csv_file["sTEC"][prn_mask & station_mask]
        jtec = csv_join["sTEC"][jprn_mask & station_mask]
        ylabel = "sTEC (TECU)"
        out_extention= "/{}-joint-sTEC_series_PRN_{}.pdf".format(station, p)
    else:
        tec = csv_file["vTEC"][prn_mask & station_mask]
        jtec = csv_join["vTEC"][jprn_mask & station_mask]


    if previous:
        TIME = np.concatenate([jtime-24.0, time])
        TEC = np.concatenate([jtec, tec])
    elif Next:
        TIME = np.concatenate([time, jtime+24.0])
        TEC = np.concatenate([tec, jtec])
    f2 = plt.figure()
    ax = f.add_subplot(n_subplots, 1, i+1)
    ax2 = f2.add_subplot(1,1,1)
    ax2.plot(TIME, TEC)
    ax2.axvline(time_hours, ls="--", c="k")
    ax2.axhline(0, ls="--", c="k")
    ax2.axvline(24., ls=":", c="b", alpha=0.5)
    ax2.set_ylim(-2.8, 2.8)
    ax2.set_xlabel("UT (hours)")
    ax2.set_ylabel(ylabel)
    f2.suptitle(station.upper()+"-PRN{}".format(p))
    f2.savefig(outdir+out_extention)
    plt.clf()
    # Now it's time to plot the collective plot
    
    ax.plot(TIME, TEC)
    ax.axvline(time_hours, ls="--", c="k")
    ax.axvline(24.0, ls=":", c="b", alpha=0.5)
    ax.set_ylim(-2.8, 2.8)
    ax.tick_params(axis="y", colors="white")
    ax.text(x_low-0.35, 0, "PRN {}".format(p), fontsize=6)
    if i < len(PRN_set)-1:
        ax.tick_params(axis="x", colors="white")
    else:
        ax.set_xlabel("Time (UT)")
    ax.set_xlim(x_low, x_high)
f.suptitle(title_label)
f.savefig(outdir+g_extention)
