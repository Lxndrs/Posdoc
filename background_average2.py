# Este programa lee el archivo creado por background_average.py, calcula la mediana del TEC, el dTEC y crea tres gráficas:
# 1. Grafica las curvas TEC de los 27 días y la mediana.
# 2. Compara la mediana con la curva TEC del día del evento
# 3. Curva de dTEC vs tiempo y en fondo de colores los correspondientes valores del indice W.

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.table import Table
from statistics import median
import pandas as pd

def julian_day(Date):
    """
    Finds the julian day of event, needed to find correctly the file name
    of the *.Std TEC data
    """
    year, month, sday = Date.split("-")
    day = float(sday)
    month_dict = {"01":[0, 0], "02":[31, 31], "03":[31+28, 31+29], "04":[31+28+31, 31+29+31], 
                  "05":[31+28+31+30, 31+29+31+30], "06":[31+28+31+30+31, 31+29+31+30+31], 
                  "07":[31+28+31+30+31+30, 31+29+31+30+31+30], "08":[31+28+31+30+31+30+31, 31+29+31+30+31+30+31], 
                  "09":[31+28+31+30+31+30+31+31, 31+29+31+30+31+30+31+31], 
                  "10":[31+28+31+30+31+30+31+31+30, 31+29+31+30+31+30+31+31+30], 
                  "11":[31+28+31+30+31+30+31+31+30+31, 31+29+31+30+31+30+31+31+30+31], 
                  "12":[31+28+31+30+31+30+31+31+30+31+30, 31+29+31+30+31+30+31+31+30+31+30]}
    bisiesto = {"2000":True, "2001":False, "2002":False, "2003":False, "2004":True, "2005":False, 
                "2006":False, "2007":False, "2008":True, "2009":False, "2010":False, "2011":False,
                "2012":True, "2013":False, "2014":False, "2015":False, "2016":True, "2017":False, 
                "2018":False, "2019":False, "2020":True, "2021":False}
    
    if bisiesto[year]:
        julian = day + month_dict[month][1]
    else:
        julian = day + month_dict[month][0]

    return str(int(julian))

# Load data from command line

parser = argparse.ArgumentParser(
	description=""" Choose a file to work""")

parser.add_argument("--station", type=str,
			default="kvtx",
			help=" Choose station")

parser.add_argument('--date', type=str, default='2000-01-01',
			help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--sets", type=str, default="1", help="Choose a set of stations")


cmd_args = parser.parse_args()
date = cmd_args.date
station = cmd_args.station
set_folder = cmd_args.sets
sns.set_style("whitegrid") # seaborn graph style
starttime = {"2019-06-22": 21+25./60.+48./3600.}


# Load previous days TEC data obtained in background_average.py
directory = "./data/{}/set{}/".format(date, set_folder)
data = Table.read(directory+"{}-Avg-prev-days-summary.tab".format(station), format="ascii")
Time = np.unique(data["time"])
ID = np.unique(data["FileID"])
TEC_median = np.zeros_like(Time)
#sTEC_m = np.zeros_like(Time)

# Load stations data for sunrise and sunset
s_data = pd.read_csv("station_data.csv")
station_mask = s_data["Site"] == station.upper()
sunrise = float(s_data["Sunrise"][station_mask])
sunset = float(s_data["Sunset"][station_mask])

# Compute medTEC and plot curves
for i,T in enumerate(Time):
    time_mask = data["time"] == T # Chhose only TEC elements that matches with the current time
    m_time = data["time"][time_mask]
    m_TEC = data["TEC"][time_mask]
#    m_sTEC = data["sigma_TEC"][time_mask] # Due to the vast number of nans in some files, this may not work
    TEC_median[i] = median(m_TEC)
#    sTEC_m[i] = max(m_sTEC) # For safety purposes, always take the maximum value of standard devation at each time
for j in ID:
    if j == len(ID) -1:
        label = "Previous days TEC curves"
    else:
        label = ""
    ID_mask = data["FileID"] ==j
    t = data["time"][ID_mask]
    tec = data["TEC"][ID_mask]
    plt.plot(t, tec, "b.", ms=0.5, label=label)
plt.plot(Time, TEC_median, "r",label= "TEC median")
plt.legend()
plt.xlabel("UT (hours)")
plt.ylabel("TEC (TECU)")
plt.title("Average TEC curves of 27 previous days compared to median, Station {}".format(station.upper()))
plt.xlim(0, 24.0)
plt.gcf().set_size_inches(8,5)
outdir = "./dtec_windex/{}/set{}/".format(date, set_folder)
plt.savefig(outdir+"{}-Avg-prev-days-summary.pdf".format(station))
plt.clf()

# Plot median against day event TEC curve

  # Load data from event 
jul = julian_day(date)
T_data = Table.read(directory+"{}{}-{}.Std".format(station.lower(), jul, date), format="ascii")
T_ime = T_data["col1"]
T_EC = T_data["col2"]
plt.plot(T_ime, T_EC, label="Event day")
plt.plot(Time, TEC_median, label="Previous days median")
plt.legend()
plt.xlabel("UT (hours)")
plt.ylabel("TEC (TECU)")
plt.title("Event day average TEC with TEC median. Station {}".format(station.upper()))
plt.xlim(0, 24.0)
plt.savefig(outdir+"{}-TEC_w_median.pdf".format(station))
plt.clf()

# Estimation of dTEC and W index graph

dTEC = np.log10(T_EC/TEC_median)
plt.plot(T_ime, dTEC)
plt.axhline(0, c="k", lw=2) # Horizontal line at dTEC=0
plt.axvline(starttime[date], ls = "--", c="k") # vertical line at fragmentation time
plt.axvline(sunrise, ls=":", c="b")
plt.axvline(sunset, ls=":", c="b")

# Shaded regions shows the W index values
# |W| < 1 --> Blue
# |W| = 2 --> Green
# |W| = 3 --> Yellow 
# |W| = 4 --> Red
plt.fill_between(T_ime, -0.046, 0.046, color="b", alpha=0.3)
plt.fill_between(T_ime, 0.046, 0.155, color="g", alpha=0.3)
plt.fill_between(T_ime, -0.046, -0.155, color="g", alpha=0.3)
plt.fill_between(T_ime, 0.155, 0.301, color="y", alpha=0.3)
plt.fill_between(T_ime, -0.155, -0.301, color="y", alpha=0.3)
plt.fill_between(T_ime, 0.301, 0.5, color="r", alpha=0.3)
plt.fill_between(T_ime, -0.301, -0.5, color="r", alpha=0.3)
plt.text(24.1, 0.35, "W=+4", color ="r")
plt.text(24.1, -0.35, "W=-4", color ="r")
plt.text(24.1, 0.22, "W=+3", color ="y")
plt.text(24.1, -0.22, "W=-3", color ="y")
plt.text(24.1, 0.1, "W=+2", color ="g")
plt.text(24.1, -0.1, "W=-2", color ="g")
plt.text(24.1, 0., r"$\mathrm{|W|<1}$", color ="w")
props = dict(boxstyle="round", facecolor="blue", alpha=0.5)
props2 = dict(boxstyle="round", facecolor="black", alpha=0.5)
plt.text(sunrise+0.1, 0.35, "Sunrise", color="w", bbox=props)
plt.text(sunset-2, 0.35, "Sunset", color="w", bbox=props)
plt.text(starttime[date]-5.5, -0.22, "Fragmentation start", color="w", bbox=props2)
plt.xlim(0, 24.0)
plt.ylim(-0.4, 0.4)
plt.xlabel("UT (hours)")
plt.ylabel(r"dTEC $\mathrm{(\log_{10}{TECU})}$")
plt.title("dTEC and W index for {} event. Station {}".format(date, station.upper()))
plt.savefig(outdir+"dTEC_Windex_{}.pdf".format(station))
