#In this program we will recover the former TEC_series.py script, since that scipt now
#plots the TEC series for several PRNs in a single graph, and I still need to see the individual graphs.

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


parser.add_argument("--previous", action="store_true", help="Plot previous day data")
parser.add_argument("--Next", action="store_true", help="Plot next day data")
parser.add_argument("--sets", type=str, help="Choose between sets of stations")
parser.add_argument("--sTEC", action="store_true", help="Use sTEC data")
cmd_args = parser.parse_args()
date = cmd_args.date
previous = cmd_args.previous
Next = cmd_args.Next
stations_set = cmd_args.sets
sTEC = cmd_args.sTEC
set_folder = "/set"+stations_set
directory = "./data/"+date+set_folder
outdir = "./TEC_series/"+date+set_folder
outlabel = ""
if previous:
    directory = "./data/{}{}/previous".format(date, set_folder)
    outlabel = "-previous"
elif Next:
    directory = "./data/{}{}/next".format(date, set_folder)
    outlabel = "-next"
if sTEC:
    infile = glob.glob(directory+"/*sTEC.csv")[0]
else:
    infile = glob.glob(directory+"/*-TEC.csv")[0]
Date, time, residual = infile.split("_")
#if previous:
#    infile = glob.glob(directory+"/previous/*.csv")[0]
#    outdir = "./TEC_series/"+date+set_folder+"/previous"
# ****************************** Read data file ************************************************

csv_file = pd.read_csv(infile)

station_array = np.unique(csv_file["Station"])
time_hours = float(time)
# fh = 130
# fl = -60
# ws = fh-fl
sns.set_style("whitegrid")

#if Next:
#    outlabel = "next"
#elif previous:
#    outlabel = "previous"
#else:
#    outlabel = ""
for s in station_array:
    s_mask = csv_file["Station"] == s
    prn_array = np.unique(csv_file["PRN"][s_mask])
    for p in prn_array:
        prn_mask = csv_file["PRN"] == p
        frequency = round(mode(np.gradient(csv_file["Time"][s_mask & prn_mask]*3600.)))
        print("Station:{}, PRN:{}, Frequency:{}".format(s, p, frequency))
        if sTEC:
            plt.plot(csv_file["Time"][prn_mask & s_mask], csv_file["sTEC"][prn_mask & s_mask])
            ylabel = "sTEC (TECU)"
            out_extention= "/{}-sTEC_series_PRN_{}{}.pdf".format(s, p, outlabel)
        else:
            plt.plot(csv_file["Time"][prn_mask & s_mask], csv_file["vTEC"][prn_mask & s_mask])
            ylabel = "vTEC (TECU)"
            out_extention = "/{}-TEC_series_PRN_{}{}.pdf".format(s, p, outlabel)
        plt.axhline(0, ls="--", c="k")
        if not Next and not previous:
            plt.axvline(time_hours, ls="--", c="k")
        plt.ylim(-2.8, 2.8)
#        plt.xticks([5, 10, 15, 20, 23], ["5:00:00", "10:00:00", "15:00:00", "20:00:00", "23:00:00"], fontsize="medium")
        plt.xlabel("UT (hours)", fontsize="large")
        plt.ylabel(ylabel, fontsize="large")
        plt.title(s.upper()+"-PRN{}".format(p), fontsize="large")
        plt.savefig(outdir+out_extention)
        plt.clf()
