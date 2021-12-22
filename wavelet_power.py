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
parser = argparse.ArgumentParser(
	description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
				   help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--frequency", type=float, default=15,
			  help="Sampling time (in seconds)")


cmd_args = parser.parse_args()
date = cmd_args.date
dt = cmd_args.frequency/3600.

directory ="./data/"+date
outdir = "./TEC_power_spectrum/"+date
pdirectory = "./data/"+date+"/previous"
infile = glob.glob(directory+"/*.csv")[0]
pinfile = glob.glob(pdirectory+"/*.csv")[0]
df = pd.read_csv(infile)
pdf = pd.read_csv(pinfile)
Date, time, residuals = infile.split("_")
start_time = float(time)
stations = np.unique(df["Station"])
levels = np.linspace(0,1,20)
norm = MidpointNormalize(midpoint=0)
for s in stations:
    station_mask = df["Station"] ==s
    PRNs = np.unique(df["PRN"][station_mask])
    prn_pairs = list(combinations(PRNs, 2))
    previous_PRNs = np.unique(pdf["PRN"][station_mask])
    previous_PRN_pairs = list(combinations(previous_PRNs, 2))
    for p,pp in zip(prn_pairs, previous_PRN_pairs):
        plt.clf()
        print("station: {}, PRNs: ({}, {})".format(s, p[0], p[1]))
        prn_mask1 = df["PRN"] ==p[0]
        prn_mask_p1 = pdf["PRN"] == pp[0]
        prn_mask2 = df["PRN"] ==p[1]
        prn_mask_p2 = pdf["PRN"] == pp[1]
        TEC_series1 = np.array(df["vTEC"][station_mask & prn_mask1])
        time_series1 = np.array(df["Time"][station_mask & prn_mask1])
        time_series2 = np.array(df["Time"][station_mask & prn_mask2])
        TEC_series2 = np.array(df["vTEC"][station_mask & prn_mask2])
        p_TEC_series1 = np.array(pdf["Time"][station_mask & prn_mask_p1])
        p_TEC_series2 = np.array(pdf["Time"][station_mask & prn_mask_p2])
        p_time_series1 = np.array(pdf["Time"][station_mask & prn_mask_p1])
        p_time_series2 = np.array(pdf["Time"][station_mask & prn_mask_p2])
#       search for intersection of time series
        c_time, c1, c2 = np.intersect1d(time_series1, time_series2, return_indices=True)
        p_ctime, pc1, pc2 = np.intersect1d(p_time_series1, p_time_series2, return_indices=True)
        w1, periods1, scales1, COI1 = wavelet(TEC_series1[c1], dt, pad=1)
        wp1, periods_p1, scales_p1, COI_p1 = wavelet(p_TEC_series1[pc1], dt, pad=1)
        w2, periods2, scales2, COI1 = wavelet(TEC_series2[c2], dt, pad=1)
        wp2, periods_p2, scales_p2, COI_p2 = wavelet(p_TEC_series2[pc2], dt, pad=1)
#       n = len(TEC_series)
#       dof = n - scales
#       signif = wave_signif(([variance]), dt=dt, scale=scales, sigtest=1, lag1=0, dof=dof)
#       sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
#       power1 = np.abs(w1)**2
#       p_power1 = np.abs(wp1)**2
#       power2 = np.abs(w2)**2
#       p_power2 = np.abs(wp2)**2 
#       sig95_norm = power/sig95
        f = plt.figure()
        ax1 = f.add_subplot(1,2,1)
        ax2 = f.add_subplot(1,2,2)
   #compute correlation function W12 = (W1 W2*)/(|W1||W2|)
        w2_conj = np.conj(w2)
        w2p_conj = np.conj(wp2)          
        w12 = (w1*w2_conj)/np.abs(w1)/np.abs(w2)
        wp12 = (wp1*w2p_conj)/np.abs(wp1)/np.abs(wp2)
     # Begin plotting
        im1=ax1.contourf(c_time, periods1, w12, levels, norm=norm, vmin=0, vmax=1)
        ax1.contour(c_time, periods1, w12, [-99, 0.95], colors=["r", "w"])
        im2 = ax2.contourf(p_ctime, periods_p1, wp12, levels, vmin=0, vmax=1) 
        ax2.contour(p_ctime, periods_p1, wp12, [-99, 0.95], colors=["r", "w"])
        ax1.set_ylim(np.min(periods1), np.max(periods1))
        ax2.set_ylim(np.min(periods_p1), np.max(periods_p1))
        ax1.invert_yaxis()
        ax2.invert_yaxis()
        ax1.fill_between(c_time, periods1[-1], COI1, facecolor=None, edgecolor="w", hatch="x", alpha=0.3)
        ax2.fill_between(p_ctime, periods_p1[-1], COI_p1, facecolor=None, edgecolor="w", hatch="x", alpha=0.3)
        ax1.plot(c_time, COI1, "w--")
        ax2.plot(p_ctime, COI_p1, "w--")
        ax1.axvline(x=start_time, c="w", ls="--")
        bar1 = f.colorbar(im1, ax=ax1)
        bar2 = f.colorbar(im2, ax=ax2)
        bar1.set_label("Correlation Degree")
        bar2.set_label("Correlation Degree")
        ax1.title.set_text("Event date")
        ax2.title.set_text("Previous day")
        ax1.set_xlabel("UT (hours)")
        ax2.set_ylabel("UT (hours)")
        ax1.set_ylabel("Period (hours)")
        ax2.set_ylabel("Period (hours)")
        plt.suptitle("TEC coherence spectrum. Date: {}. Station: {}. PRNs: {}, {}".format(date, s, p[0], p[1]))
        #plt.title("{}-PRN{}. Date: {}".format(s, p, date))
        #plt.xlabel("UT (hours)")
        #plt.ylabel("Period (hours)")
        f.set_size_inches(18, 8)
        plt.savefig(outdir+"/{}-PRNs{}_{}-contour.pdf".format(s, p[0],p[1]))
