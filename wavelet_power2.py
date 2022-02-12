#We used a github repository to analyze the wavelets spectrum base on Torrence & Compo 1998
# Available at https://github.com/chris-torrence/wavelets
# This is a clone of wavelet_power.py with the exception we will get the wavelet spectrum for a single station and a single PRN
# (No coherence spectrum) with a contour of the 95% confidence interval
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
parser.add_argument("--station", type=str, help="Choose station")
parser.add_argument("--PRN", type=int, help="Choose PRN")
cmd_args = parser.parse_args()
date = cmd_args.date
dt = cmd_args.frequency/3600.
stations_set = cmd_args.sets
set_folder = "/set"+stations_set
station = cmd_args.station
PRN = cmd_args.PRN
directory ="./data/"+date+set_folder
outdir = "./paper/figures"
pdirectory = "./data/"+date+set_folder+"/previous"
infile = glob.glob(directory+"/*.csv")[0]
pinfile = glob.glob(pdirectory+"/*.csv")[0]
df = pd.read_csv(infile)
pdf = pd.read_csv(pinfile)
Date, time, residuals = infile.split("_") 
start_time = float(time)
levels = 20
norm = MidpointNormalize(midpoint=0)

fh = 130
fl = -60
ws = fh-fl

prn_mask = df["PRN"] == PRN
p_prn_mask = pdf["PRN"] == PRN
s_mask = df["Station"] == station
ps_mask = pdf["Station"] == station
TEC_series = np.array(df["vTEC"][prn_mask & s_mask])
p_TEC_series = np.array(pdf["vTEC"][p_prn_mask & ps_mask])
time_series = df["Time"][prn_mask & s_mask]
p_time_series = pdf["Time"][p_prn_mask & ps_mask]
w, periods, scales, COI = wavelet(TEC_series, dt, pad=1)
wp, periods_p, scales_p, COI_p = wavelet(p_TEC_series, dt, pad=1)
# Begin plotting
f = plt.figure()
ax1 = f.add_subplot(1,2,1)
ax2 = f.add_subplot(1,2,2)
power = (np.abs(w))**2
p_power = (np.abs(wp))**2
max_scale = max(np.max(power), np.max(p_power))
im1 = ax1.contourf(time_series, periods, power, levels, vmin=0, vmax=max_scale)
#ax1.contour(time_series, periods, (np.abs(w))**2, [-99, 0.95], colors=["r", "w"])
im2 = ax2.contourf(p_time_series, periods_p, (np.abs(wp))**2, levels, vmin=0, vmax=max_scale)
#ax2.contour(p_time_series, periods_p, (np.abs(wp))**2, [-99, 0.95], colors=["r", "w"])
ax1.set_ylim(np.min(periods), np.max(periods))
ax2.set_ylim(np.min(periods_p), np.max(periods_p))    
ax1.invert_yaxis()
ax2.invert_yaxis()
ax1.fill_between(time_series, periods[-1], COI, facecolor=None, edgecolor="w", hatch="x", alpha=0.3)
ax2.fill_between(p_time_series, periods_p[-1], COI_p, facecolor=None, edgecolor="w", hatch="x", alpha=0.3)
ax1.plot(time_series, COI, "w--")
ax2.plot(p_time_series, COI_p, "w--")
ax1.axvline(x=start_time, c="w", ls="--")
bar1 = f.colorbar(im1, ax=ax1)
bar2 = f.colorbar(im2, ax=ax2)
bar1.set_label("Wavelet amplitude")
bar2.set_label("Wavelet amplitude")
ax1.set_xlabel("UT (hours)")
ax1.set_ylabel("Period (hours)")
ax1.title.set_text("Event date")  
ax2.title.set_text("Previous day")
plt.suptitle("TEC power spectrum. Date: {}. Station: {}. PRN: {}".format(date, station, PRN))
f.set_size_inches(18, 8)
plt.savefig(outdir+"/{}-PRN{}-contour.pdf".format(station, PRN))
