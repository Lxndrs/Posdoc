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
parser.add_argument("--sTEC", action="store_true", help="Choose sTEC data")
parser.add_argument("--previous", action="store_true", help="Use previous day data")
parser.add_argument("--Next", action="store_true", help="Use next day data")
cmd_args = parser.parse_args()
date = cmd_args.date
dt = cmd_args.frequency/3600.
stations_set = cmd_args.sets
set_folder = "/set"+stations_set
station = cmd_args.station.lower()
PRN = cmd_args.PRN
previous = cmd_args.previous
Next = cmd_args.Next
directory ="./data/"+date+set_folder
outlabel=""
if previous:
    directory = directory +"/previous"
    outlabel = "-previous"
elif Next:
    directory = directory + "/Next"
    outlabel = "-next"
outdir = "./power_spec/{}/set{}/".format(date,stations_set)
#pdirectory = "./data/"+date+set_folder+"/previous"
infile = glob.glob(directory+"/*.csv")[0]
#pinfile = glob.glob(pdirectory+"/*.csv")[0]
df = pd.read_csv(infile)
#pdf = pd.read_csv(pinfile)
Date, time, residuals = infile.split("_") 
#start_time = float(time)
levels = 20
norm = MidpointNormalize(midpoint=0)
sTEC = cmd_args.sTEC
fh = 130
fl = -60
ws = fh-fl

prn_mask = df["PRN"] == PRN
#p_prn_mask = df["PRN"] == PRN
s_mask = df["Station"] == station
#ps_mask = pdf["Station"] == station
if sTEC:
    TEC_series = np.array(df["sTEC"][prn_mask & s_mask])
#    p_TEC_series = np.array(pdf["sTEC"][p_prn_mask & ps_mask])
else:
    TEC_series = np.array(df["vTEC"][prn_mask & s_mask])
#    p_TEC_series = np.array(pdf["vTEC"][p_prn_mask & ps_mask])
time_series = df["Time"][prn_mask & s_mask]
#p_time_series = pdf["Time"][p_prn_mask & ps_mask]
w, periods, scales, COI = wavelet(TEC_series, dt, pad=1)
#wp, periods_p, scales_p, COI_p = wavelet(p_TEC_series, dt, pad=1)
# Begin plotting
f = plt.figure()
ax1 = f.add_subplot(1,1,1)
#ax2 = f.add_subplot(1,2,2)
power = (np.abs(w))**2
#p_power = (np.abs(wp))**2
#max_scale = max(np.max(power), np.max(p_power))
freq = 1e3/(3600*periods) 
im1 = ax1.contourf(time_series, freq, power, levels)
#print(im1.levels)
#im2 = ax2.contourf(p_time_series, periods_p, (np.abs(wp))**2, levels, vmin=0, vmax=max_scale)
var = TEC_series - np.mean(TEC_series)
variance = np.std(var, ddof=1)**2
signif = wave_signif(([variance]), dt=dt, sigtest=0, scale=scales, lag1=0, mother="MORLET")
n = len(TEC_series)
sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
sig_95 = power/sig95
ax1.contour(time_series, freq, sig_95, [-99, 1], colors=["r", "w"])
ax1.set_ylim(np.min(freq), 2)
#ax2.set_ylim(np.min(periods_p), np.max(periods_p))    
#ax1.invert_yaxis()
#ax2.invert_yaxis()
ax1.fill_between(time_series, freq[-1], 1e3/(3600*COI), facecolor=None, edgecolor="w", hatch="x", alpha=0.3)
#ax2.fill_between(p_time_series, periods_p[-1], COI_p, facecolor=None, edgecolor="w", hatch="x", alpha=0.3)
ax1.plot(time_series, 1e3/(3600*COI), "w--")
#ax2.plot(p_time_series, COI_p, "w--")
#ax1.axvline(x=start_time, c="w", ls="--")
bar1 = f.colorbar(im1, ax=ax1)
#bar2 = f.colorbar(im2, ax=ax2)
bar1.set_label("Wavelet amplitude")
#bar2.set_label("Wavelet amplitude")
ax1.set_xlabel("UT (hours)")
ax1.set_ylabel("Frequency (mHz)")
ax1.title.set_text("Date: {}. Station: {}. PRN: {}".format(date, station.upper(), PRN))  
#ax2.title.set_text("Previous day")
ax1.set_xlim(min(time_series), max(time_series))
dat0=im1.allsegs[-1]
for seg in dat0:
    plt.text(min(time_series)+1, 1.9, "time={:.3f}".format(np.mean(seg[:,0])), color="white")
    plt.text(min(time_series)+1, 1.8, "frequency ={:.3f}".format(np.mean(seg[:,1])), color="white")
plt.suptitle("TEC power spectrum")
f.set_size_inches(10, 10)
if sTEC:
    plt.savefig(outdir+"/{}-PRN{}-sTEC-contour{}.pdf".format(station, PRN, outlabel))
else:
    plt.savefig(outdir+"/{}-PRN{}-contour.pdf{}".format(station, PRN, outlabel))
