#We used a github repository to analyze the wavelets spectrum base on Torrence & Compo 1998
# Available at https://github.com/chris-torrence/wavelets
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

parser.add_argument("--sets", type=str, help="Choose between sets of stations")
cmd_args = parser.parse_args()
date = cmd_args.date
dt = cmd_args.frequency/3600.
stations_set = cmd_args.sets
set_folder = "/set"+stations_set
directory ="./data/"+date+set_folder
outdir = "./TEC_power_spectrum/"+date+set_folder
pdirectory = "./data/"+date+set_folder+"/previous"
infile = glob.glob(directory+"/*.csv")[0]
pinfile = glob.glob(pdirectory+"/*.csv")[0]
df = pd.read_csv(infile)
pdf = pd.read_csv(pinfile)
Date, time, residuals = infile.split("_") 
start_time = float(time)
stations = np.unique(df["Station"])
#window_size = np.unique(df["WSize"])
levels = np.linspace(0,1,20)
norm = MidpointNormalize(midpoint=0)
station_pairs = list(combinations(stations, 2))
#for ws in window_size: # This should be working as follows: PRNs --> Stations (Pair stations with same PRNs)
#    win_mask = df["WSize"] == ws
#    p_win_mask = pdf["WSize"] == ws
fh = 130
fl = -60
ws = fh-fl
PRNs = np.unique(df["PRN"])
for p in PRNs:
      prn_mask = df["PRN"] == p
      p_prn_mask = pdf["PRN"] == p
      for s1, s2 in station_pairs:
          plt.clf()
          print("Stations: {} {}, PRN: {}".format(s1, s2, p))
          s_mask1 = df["Station"] == s1
          ps_mask1 = pdf["Station"] == s1
          TEC_series1 = np.array(df["vTEC"][prn_mask & s_mask1])
          p_TEC_series1 = np.array(pdf["vTEC"][p_prn_mask & ps_mask1])
          s_mask2 = df["Station"] == s2
          ps_mask2 = pdf["Station"] == s2
          TEC_series2 = np.array(df["vTEC"][prn_mask & s_mask2])
          p_TEC_series2 = np.array(pdf["vTEC"][prn_mask & ps_mask2])
          time_series1 = df["Time"][prn_mask & s_mask1]
          time_series2 = df["Time"][prn_mask & s_mask2]
          p_time_series1 = pdf["Time"][prn_mask & ps_mask1]
          p_time_series2 = pdf["Time"][prn_mask & ps_mask2]
          # Ignore combinations where at least one array is empty
          if (len(time_series1)==0) | (len(time_series2)==0):
              print("Times series empty")
              continue
          if (len(p_time_series1)==0) | (len(p_time_series2)==0):
              print("Times series empty for previous day")
              continue
          #search for intersection of time series                                                    
          c_time, c1, c2 = np.intersect1d(time_series1, time_series2, return_indices=True)       
          p_ctime, pc1, pc2 = np.intersect1d(p_time_series1, p_time_series2, return_indices=True)
          # compute wavelet transforms
          try:
              w1, periods1, scales1, COI1 = wavelet(TEC_series1[c1], dt, pad=1)
              w2, periods2, scales2, COI2 = wavelet(TEC_series2[c2], dt, pad=1)
              wp1, periodsp1, scalesp1, COIp1 = wavelet(p_TEC_series1[pc1], dt, pad=1)
              wp2, periodsp2, scalesp2, COIp2 = wavelet(p_TEC_series2[pc2], dt, pad=1)
          except ValueError:
              print("No match between PRNs")
              continue
          w2_conj = np.conj(w2)
          wp2_conj = np.conj(wp2)
          #compute correlation function W12 = (W1 W2*)**2/(|W1||W2|)**2
          w12 = (np.abs(w1*w2_conj)/np.abs(w1)/np.abs(w2))
          wp12 = (np.abs(wp1*wp2_conj)/np.abs(wp1)/np.abs(wp2))
          # Begin plotting
          f = plt.figure()
          ax1 = f.add_subplot(1,2,1)
          ax2 = f.add_subplot(1,2,2)
          im1 = ax1.contourf(c_time, periods1, w12, levels, vmin=0, vmax=1)
          ax1.contour(c_time, periods1, w12, [-99, 0.95], colors=["r", "w"])
          im2 = ax2.contourf(p_ctime, periodsp1, wp12, levels, vmin=0, vmax=1)
          ax2.contour(p_ctime, periodsp1, wp12, [-99, 0.95], colors=["r", "w"])
          ax1.set_ylim(np.min(periods1), np.max(periods1))
          ax2.set_ylim(np.min(periodsp1), np.max(periodsp1))    
          ax1.invert_yaxis()
          ax2.invert_yaxis()
          ax1.fill_between(c_time, periods1[-1], COI1, facecolor=None, edgecolor="w", hatch="x", alpha=0.3)
          ax1.plot(c_time, COI1, "w--")
          ax2.plot(p_ctime, COIp1, "w--")
          ax1.axvline(x=start_time, c="w", ls="--")
          bar1 = f.colorbar(im1, ax=ax1)
          bar2 = f.colorbar(im2, ax=ax2)
          bar1.set_label("Correlation Degree")
          bar2.set_label("Correlation Degree")
          ax1.set_xlabel("UT (hours)")
          ax1.set_ylabel("Period (hours)")
          ax1.title.set_text("Event date")  
          ax2.title.set_text("Previous day")
          plt.suptitle("TEC coherence spectrum. Date: {}. Stations: {},{}. PRN: {}".format(date, s1, s2, p))
          f.set_size_inches(18, 8)
          plt.savefig(outdir+"/{} {}-PRN{}-contour.pdf".format(s1, s2, p))
