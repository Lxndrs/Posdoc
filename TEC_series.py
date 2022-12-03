import pandas as pd
import argparse 
import numpy as np
import matplotlib.pyplot as plt
import glob
from statistics import mode
import seaborn as sns
parser = argparse.ArgumentParser(description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
				      help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--station", type=str, default="cn05", help="Choose station")
parser.add_argument("--PRNs", type=str, default="1", help="Choose PRNs to graph (separate with commas)")
parser.add_argument("--previous", action="store_true", help="Choose previous day data")
parser.add_argument("--Next", action="store_true", help="Choose next day data")
parser.add_argument("--sets", type=str, default="1", help="Choose between sets of stations")
parser.add_argument("--allp", action="store_true", help="Use all PRNs in graph")
parser.add_argument("--sTEC", action="store_true", help="Use sTEC data")
cmd_args = parser.parse_args()
date = cmd_args.date
previous = cmd_args.previous
Next = cmd_args.Next
stations_set = cmd_args.sets
set_folder = "/set"+stations_set
directory = "./data/"+date+set_folder
outdir = "./TEC_series/"+date+set_folder
prn_str = cmd_args.PRNs
s = cmd_args.station.lower()
allp = cmd_args.allp
sTEC = cmd_args.sTEC
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
# ****************************** Read data file ************************************************

csv_file = pd.read_csv(infile)
time_hours = float(time)
x_low = max(time_hours - 1.1, 0.0)
x_high = min(time_hours + 3.35, 24.0)
sns.set_style("whitegrid")
s_mask = csv_file["Station"] == s
if allp:
    prn_array = np.unique(csv_file["PRN"][s_mask])
else:
    prn_array_str = prn_str.split(",")
    prn_array = [int(i) for i in prn_array_str]
n_subplots = len(prn_array)
fig = plt.figure()
big_ax = fig.add_subplot(1,1,1)
# Turn off axis lines and ticks of the big subplot
big_ax.spines['top'].set_color('none')
big_ax.spines['bottom'].set_color('none')
big_ax.spines['left'].set_color('none') 
big_ax.spines['right'].set_color('none')
big_ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
for i, p in enumerate(prn_array):
    ax = fig.add_subplot(n_subplots, 1, i+1)
    prn_mask = csv_file["PRN"] == p
#    frequency = round(mode(np.gradient(csv_file["Time"][s_mask & prn_mask]*3600.)))
    print("Station:{}, PRN:{}".format(s, p))
    if sTEC:
        ax.plot(csv_file["Time"][prn_mask & s_mask], csv_file["sTEC"][prn_mask & s_mask])
        out_extention = "/{}-sTEC_series{}.pdf".format(s, outlabel)
    else:  
        ax.plot(csv_file["Time"][prn_mask & s_mask], csv_file["vTEC"][prn_mask & s_mask])
        out_extention = "/{}-TEC_series{}.pdf".format(s, outlabel)
    if previous:
        x_low = 22.0
        x_high = 24.1
    elif Next:
        x_low = 0.0
        x_high = 8.0
    else:
        ax.axvline(time_hours, ls="--", c="k")
    plt.ylim(-2.8, 2.8)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(axis="x", labelsize=8)
    #ax.set_ylabel("sTEC (TECU)", fontsize="x-small")
    ax.text(x_low, 2.82, "PRN {}".format(p), fontsize=6)

    if i < len(prn_array)-1:
        ax.tick_params(axis="x", colors="white")
    ax.set_xlim(x_low, x_high)
big_ax.set_xlabel("Time (UT)")
big_ax.grid(False)
big_ax.set_ylabel("sTEC (TECU)")
if sTEC:
    fig.suptitle("Detrended slant TEC series for individual satellites. Station {}".format(s.upper()))
else:
    fig.suptitle("Detrended vTEC series for individual satellites. Station {}".format(s.upper()))
plt.savefig(outdir+out_extention)
