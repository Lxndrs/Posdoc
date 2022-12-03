#We used a github repository to analyze the wavelets spectrum base on Torrence & Compo 1998
# Available at https://github.com/chris-torrence/wavelets
#This program is intended to calculate coherence spectrums between two time series, we must select them by
# giving the program the corresponding station and satellite which collected the data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
sys.path.insert(0, "../wavelets/wave_python/") 
from waveletFunctions import wavelet, wave_signif
import matplotlib.colors as mcolors
from midpoint import MidpointNormalize
import glob
from itertools import combinations
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser(
	  description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
				     help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--frequency", type=float, default=15,
			    help="Sampling time (in seconds)")

parser.add_argument("--sets1", type=str, default="1", help="Choose between sets of stations (first data set)")
parser.add_argument("--sets2", type=str, default="1", help="Choose between sets of stations (second data set)")
parser.add_argument("--station1", type=str, default="bara", help="Choose first station")
parser.add_argument("--station2", type=str, default="bara", help="Choose second station")
parser.add_argument("--PRN1", type=int, default=1, help="Select the first satellite PRN")
parser.add_argument("--PRN2", type=int, default=1, help="Select the second satellite PRN")
parser.add_argument("--previous", action="store_true", help="Use previous day data")
parser.add_argument("--Next", action="store_true", help="Use next day data")
parser.add_argument("--sTEC", action="store_true", help="Use sTEC data")
cmd_args = parser.parse_args()
date = cmd_args.date
dt = cmd_args.frequency/3600.
stations_set1 = cmd_args.sets1
set_folder1 = "/set"+stations_set1
stations_set2 = cmd_args.sets2
set_folder2 = "/set"+stations_set2
directory1 ="./data/"+date+set_folder1
directory2 ="./data/"+date+set_folder2
outdir = "./TEC_power_spectrum/"+date+"/outputs"
outlabel=""
previous = cmd_args.previous
Next = cmd_args.Next
if previous:
    directory1 = directory1+"/previous"
    directory2 = directory2+"/previous"
    outlabel = "-previous"
elif Next:
    directory1 = directory1+"/next"
    directory2 = directory2+"/next"
    outlabel = "-next"
s1 = cmd_args.station1
s2 = cmd_args.station2
p1 = cmd_args.PRN1
p2 = cmd_args.PRN2
infile1 = glob.glob(directory1+"/*.csv")[0]
infile2 = glob.glob(directory2+"/*.csv")[0]
sTEC = cmd_args.sTEC
#pinfile = glob.glob(pdirectory+"/*.csv")[0]
df1 = pd.read_csv(infile1)
df2 = pd.read_csv(infile2)
#Date, time, residuals = infile.split("_") 
#start_time = float(time)
#stations = np.unique(df["Station"])
#window_size = np.unique(df["WSize"])
levels = np.linspace(0,1,20)
norm = MidpointNormalize(midpoint=0)
#station_pairs = list(combinations(stations, 2))
#for ws in window_size: # This should be working as follows: PRNs --> Stations (Pair stations with same PRNs)
#    win_mask = df["WSize"] == ws
#    p_win_mask = pdf["WSize"] == ws
#fh = 130
#fl = -60
#ws = fh-fl
#PRNs = np.unique(df["PRN"])
#for p in PRNs:
prn_mask1 = df1["PRN"] == p1
prn_mask2 = df2["PRN"] == p2
#      for s1, s2 in station_pairs:
#          plt.clf()
#          print("Stations: {} {}, PRN: {}".format(s1, s2, p))
s_mask1 = df1["Station"] == s1
s_mask2 = df2["Station"] == s2
if sTEC:
    TEC_series1 = np.array(df1["sTEC"][prn_mask1 & s_mask1])
    TEC_series2 = np.array(df2["sTEC"][prn_mask2 & s_mask2])
else:
    TEC_series1 = np.array(df1["vTEC"][prn_mask1 & s_mask1])
    TEC_series2 = np.array(df2["vTEC"][prn_mask2 & s_mask2])
 #         s_mask2 = df["Station"] == s2
#          ps_mask2 = pdf["Station"] == s2
#          TEC_series2 = np.array(df["vTEC"][prn_mask & s_mask2])
#          p_TEC_series2 = np.array(pdf["vTEC"][prn_mask & ps_mask2])
time_series1 = df1["Time"][prn_mask1 & s_mask1]
time_series2 = df2["Time"][prn_mask2 & s_mask2]
#          p_time_series1 = pdf["Time"][prn_mask & ps_mask1]
#          p_time_series2 = pdf["Time"][prn_mask & ps_mask2]
          # Ignore combinations where at least one array is empty
#          if (len(time_series1)==0) | (len(time_series2)==0):
#              print("Times series empty")
#              continue
#          if (len(p_time_series1)==0) | (len(p_time_series2)==0):
#              print("Times series empty for previous day")
#              continue
#search for intersection of time series                                                    
c_time, c1, c2 = np.intersect1d(time_series1, time_series2, return_indices=True)       
#p_ctime, pc1, pc2 = np.intersect1d(p_time_series1, p_time_series2, return_indices=True)
# compute wavelet transforms
try:
    w1, periods1, scales1, COI1 = wavelet(TEC_series1[c1], dt, pad=1)
    w2, periods2, scales2, COI2 = wavelet(TEC_series2[c2], dt, pad=1)
 #             wp1, periodsp1, scalesp1, COIp1 = wavelet(p_TEC_series1[pc1], dt, pad=1)
 #             wp2, periodsp2, scalesp2, COIp2 = wavelet(p_TEC_series2[pc2], dt, pad=1)
except ValueError:
    print("No match between PRNs")
    sys.exit()
w2_conj = np.conj(w2)
#wp2_conj = np.conj(wp2)
#compute correlation function W12 = (W1 W2*)**2/(|W1||W2|)**2
w12 = (np.abs(w1*w2_conj)/np.abs(w1)/np.abs(w2))
#wp12 = (np.abs(wp1*wp2_conj)/np.abs(wp1)/np.abs(wp2))
# Begin plotting
#          f = plt.figure()
#          ax1 = f.add_subplot(1,2,1)
#          ax2 = f.add_subplot(1,2,2)
im1 = plt.contourf(c_time, periods1, w12, levels, vmin=0, vmax=1)
plt.contour(c_time, periods1, w12, [-99, 0.95], colors=["r", "w"])
#          im2 = ax2.contourf(p_ctime, periodsp1, wp12, levels, vmin=0, vmax=1)
#          ax2.contour(p_ctime, periodsp1, wp12, [-99, 0.95], colors=["r", "w"])
plt.ylim(np.min(periods1), np.max(periods1))
#          ax2.set_ylim(np.min(periodsp1), np.max(periodsp1))    
plt.gca().invert_yaxis()
#          ax2.invert_yaxis()
plt.fill_between(c_time, periods1[-1], COI1, facecolor=None, edgecolor="w", hatch="x", alpha=0.3)
plt.plot(c_time, COI1, "w--")
#ax2.plot(p_ctime, COIp1, "w--")
#plt.axvline(x=start_time, c="w", ls="--")
bar1 = plt.colorbar(im1)
#bar2 = plt.colorbar(im2, ax=ax2)
bar1.set_label("Correlation Degree")
#bar2.set_label("Correlation Degree")
plt.xlabel("UT (hours)")
plt.ylabel("Period (hours)")
#plt.title("Event date")  
#ax2.title.set_text("Previous day")
plt.suptitle("TEC coherence spectrum. Date: {}. Stations: {},{}. PRNs: {}, {}".format(date, s1, s2, p1, p2))
#          f.set_size_inches(18, 8)
plt.savefig(outdir+"/{} {}-PRNs{}-{}-contour{}.pdf".format(s1, s2, p1, p2, outlabel))
