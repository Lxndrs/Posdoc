# Este programa lee el archivo creado por background_average.py, calcula la mediana del TEC, el dTEC y crea las siguiente grafica:
# 1. Curva de dTEC vs tiempo y en fondo de colores los correspondientes valores del indice W para todas las estaciones que NO hayan
# detectado TIDs y calcula la mediana de todas estas curvas.
# 2. Resta del dTEC de las estaciones que si realizaron deteccion con la mediana

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.table import Table
from statistics import median, quantiles
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

#parser.add_argument("--station", type=str,
#			default="kvtx",
#			help=" Choose station")

parser.add_argument('--date', type=str, default='2000-01-01',
			help='Choose date. Format: yyyy-mm-dd')

#parser.add_argument("--sets", type=str, default="1", help="Choose a set of stations")


cmd_args = parser.parse_args()
date = cmd_args.date
#station = cmd_args.station
#set_folder = cmd_args.sets
sns.set_style("whitegrid") # seaborn graph style
starttime = {"2019-06-22": 21+25./60.+48./3600.}
outdir = "./dtec_windex/{}/".format(date)
jul = julian_day(date)
# Enumerate stations list with non detections

n_stations = ["BARA", "CN05", "CN27", "CRLR", "CRSE", "JME2", "LVEG", #"RDAZ", 
"RDF2", "RDHI", "RDLT", "RDMA", "RDMC", "RDMS", "RDNE", "RDSD", "RDSF", "RDSJ", "SPED", 
"SROD", "TGDR", "AIRS", "CN00", "GERD", "NWBL", "OLVN", "RCHY", "RDON"]#, "TRNT"]
n_sets = {"BARA":1, "CN05":1, "CN27":1, "CRLR":1, "CRSE":1, "JME2":1, "LVEG":1, #"RDAZ":1, 
"RDF2":1, "RDHI":1, "RDLT":1, "RDMA":1, "RDMC":1, "RDMS":1, "RDNE":1, "RDSD":1, "RDSF":1, 
"RDSJ":1, "SPED":1, "SROD":1, "TGDR":1, "AIRS":2, "CN00":2, "GERD":2, "NWBL":2, "OLVN":2, 
"RCHY":2, "RDON":2}#, "TRNT":2}
#RDAZ and TRNT stations somehow make the program crash
# Enumerate stations with detections

stations = ["CN04", "GRE1", "CN19", "CN40", "TTSF", "TTUW", "BOAV", "KOUG", "KOUR"] 
sets= {"CN04":3, "GRE1":3, "CN19":4, "CN40":4, "TTSF":4, "TTUW":4, "BOAV":5, "KOUG":5, "KOUR":5}
# Set output table for computing median dTEC

outtab = Table(names=("Time", "dTEC"))

# Load previous days TEC data obtained in background_average.py and iterate in n_stations
for ns in n_stations:
#    print("Station: {}".format(ns))
    directory = "./data/{}/set{}/".format(date, n_sets[ns])
    data = Table.read(directory+"{}-Avg-prev-days-summary.tab".format(ns.lower()), format="ascii")
    Time = np.unique(data["time"])
    ID = np.unique(data["FileID"])
    TEC_median = np.zeros_like(Time)

# Load data from event
    T_data = Table.read(directory+"{}{}-{}.Std".format(ns.lower(), jul, date), format="ascii")
    T_ime = T_data["col1"]
    T_EC = T_data["col2"]
    dTEC = np.zeros_like(T_ime)

    
# Compute medTEC 
    for i,T in enumerate(Time):
        time_mask = data["time"] == T # Choose only TEC elements that matches with the current time
        m_time = data["time"][time_mask]
        m_TEC = data["TEC"][time_mask]
        TEC_median[i] = median(m_TEC)

# Estimate  dTECs
        dTEC[i] = np.log10(T_EC[i]/TEC_median[i])
        outtab.add_row([T, dTEC[i]])
#    plt.plot(T_ime, dTEC, "b", lw=0.5)


#    plt.plot(T_ime, d_TEC)

time = np.unique(outtab["Time"])
DTEC = np.unique(outtab["dTEC"])
dTEC_median = []
#dTEC_Q1 = []
#dTEC_Q3 = []
for t, dt in zip(time, DTEC):
    t_mask = outtab["Time"] == t
    mdtec = outtab["dTEC"][t_mask]
    dtec_median = median(mdtec)
    dTEC_median.append(dtec_median)
#    q1, q2, q3 = quantiles(mdtec)
#    dTEC_Q1.append(q1)
#    dTEC_Q3.append(q3)

# Load stations data for sunrise and sunset
st_data = pd.read_csv("station_data.csv")

# Properties of text boxes
props = dict(boxstyle="round", facecolor="blue", alpha=0.5)
props2 = dict(boxstyle="round", facecolor="black", alpha=0.5)

# Plot TEC median 
plt.plot(time, dTEC_median, "r")
plt.xlim(0, 24.0)
plt.ylim(-0.4, 0.4)
plt.xticks([5, 10, 15, 20, 23], ["5:00:00", "10:00:00", "15:00:00", "20:00:00", "23:00:00"])
plt.xlabel("UT (hours)", fontsize="large")
plt.ylabel("dTEC median", fontsize="large")
plt.title("Median dTEC (of 'non detecting' stations)", fontsize="large")
plt.gcf().set_size_inches(8,5)
plt.savefig(outdir+"/ddtec/ddTEC-median.pdf")
plt.clf()
# Estimate and plot dTEC median and station dTEC
for s in stations:
    s_data = Table.read(outdir+"set{}/dTEC_table_{}.csv".format(sets[s], s.lower()), format="csv")
    T_ime = s_data["Time"]
    d_TEC = s_data["dTEC"]
    station_mask = st_data["Site"] == s.upper()
    sunrise = float(st_data["Sunrise"][station_mask])
    sunset = float(st_data["Sunset"][station_mask])
    plt.plot(time, d_TEC - dTEC_median)
    plt.axhline(0, c="k", lw=2) # Horizontal line at dTEC=0
    plt.axvline(starttime[date], ls = "--", c="k") # vertical line at fragmentation time
    plt.axvline(sunrise, ls=":", c="b")
    plt.axvline(sunset, ls=":", c="b")
    plt.text(sunrise+0.1, 0.35, "Sunrise", color="w", bbox=props, fontsize="medium")
    plt.text(sunset-2, 0.35, "Sunset", color="w", bbox=props, fontsize="medium")
    plt.text(starttime[date]-5.5, -0.22, "Fragmentation start", color="w", bbox=props2, fontsize="medium")
    plt.xlim(0, 24.0)
    plt.ylim(-0.4, 0.4)
    plt.xticks([5, 10, 15, 20, 23], ["5:00:00", "10:00:00", "15:00:00", "20:00:00", "23:00:00"])
    plt.xlabel("UT (hours)", fontsize="large")
    plt.ylabel("ddTEC", fontsize="large")
    plt.title("dTEC of station {} minus median dTEC (of 'non detecting' stations).".format(s), fontsize="large")
    plt.gcf().set_size_inches(8,5)
    plt.savefig(outdir+"/ddtec/ddTEC-{}.pdf".format(s))
    plt.clf()
#print(len(time), len(dTEC_Q1), len(dTEC_Q3))
#IQR = np.array(dTEC_Q3) - np.array(dTEC_Q1)



# Load stations data for sunrise and sunset
#s_data = pd.read_csv("station_data.csv")
#station_mask = s_data["Site"] == station.upper()
#sunrise = float(s_data["Sunrise"][station_mask])
#sunset = float(s_data["Sunset"][station_mask])

# Load data from event 
#
#
#
#

# Estimation of dTEC and W index graph

#
#plt.plot(T_ime, dTEC)
#plt.axhline(0, c="k", lw=2) # Horizontal line at dTEC=0
#plt.axvline(starttime[date], ls = "--", c="k") # vertical line at fragmentation time
#plt.axvline(sunrise, ls=":", c="b")
#plt.axvline(sunset, ls=":", c="b")

# Shaded regions shows the W index values
# |W| < 1 --> Blue
# |W| = 2 --> Green
# |W| = 3 --> Yellow 
# |W| = 4 --> Red
#plt.fill_between(T_ime, -0.046, 0.046, color="b", alpha=0.3)
#plt.fill_between(T_ime, 0.046, 0.155, color="g", alpha=0.3)
#plt.fill_between(T_ime, -0.046, -0.155, color="g", alpha=0.3)
#plt.fill_between(T_ime, 0.155, 0.301, color="y", alpha=0.3)
#plt.fill_between(T_ime, -0.155, -0.301, color="y", alpha=0.3)
#plt.fill_between(T_ime, 0.301, 0.5, color="r", alpha=0.3)
#plt.fill_between(T_ime, -0.301, -0.5, color="r", alpha=0.3)
#plt.text(24.1, 0.35, "W=+4", color ="r")
#plt.text(24.1, -0.35, "W=-4", color ="r")
#plt.text(24.1, 0.22, "W=+3", color ="y")
#plt.text(24.1, -0.22, "W=-3", color ="y")
#plt.text(24.1, 0.1, "W=+2", color ="g")
#plt.text(24.1, -0.1, "W=-2", color ="g")
#plt.text(24.1, 0., r"$\mathrm{|W|<1}$", color ="b")
#props = dict(boxstyle="round", facecolor="blue", alpha=0.5)
#props2 = dict(boxstyle="round", facecolor="black", alpha=0.5)
#plt.text(sunrise+0.1, 0.35, "Sunrise", color="w", bbox=props, fontsize="medium")
#plt.text(sunset-2, 0.35, "Sunset", color="w", bbox=props, fontsize="medium")
#plt.text(starttime[date]-5.5, -0.22, "Fragmentation start", color="w", bbox=props2, fontsize="medium")
#plt.xticks([5, 10, 15, 20, 23], ["5:00:00", "10:00:00", "15:00:00", "20:00:00", "23:00:00"], fontsize="medium")
#plt.xlabel("UT (hours)", fontsize="large")
#plt.ylabel(r"dTEC $\mathrm{(\log_{10}{TECU})}$", fontsize="large")
#plt.title("dTEC and W index for {} event. Station {}".format(date, station.upper()), fontsize="large")
#plt.savefig(outdir+"dTEC_Windex_all{}.pdf".format(station))
