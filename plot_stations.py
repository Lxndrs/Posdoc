# Mexico map plotter
# The main idea of this program was taken from 
# https://towardsdatascience.com/mapping-with-matplotlib-pandas-geopandas-and-basemap-in-python-d11b57ab5dac
# By Ashwani Dhankhar 
# And the shape file for Mexico from CONABIO
# http://www.conabio.gob.mx/informacion/metadata/gis/destdv250k_2gw.xml?_xsl=/db/meadata/xsl/fgdc_html.xsl&_indent=no

# The goal of this program was settled to plot event trajectory from GLM data and the stations that were used to get the RINEX data 

import seaborn as sns
import numpy as np
import pandas as pd
import shapefile as shp
import matplotlib.pyplot as plt
from plotfullmap import plot_map
import argparse
from astropy.table import Table
import glob
import argparse

# Load desired date from command line

parser = argparse.ArgumentParser(
  description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
				 help='Choose date. Format: yyyy-mm-dd')



cmd_args = parser.parse_args()
date = cmd_args.date


directory = "./data/"+date+"/GLM/"


# set figure style
sns.set_style("whitegrid") 
sns.mpl.rc("figure", figsize=(10,6))

# Read shape file of Mexico map
sf = shp.Reader("map.shp")
plot_map(sf)

# Read stations positions table

stations_pos = Table.read("station_data.tab", format="ascii")

# Load and plot event position
load_meteor_pos = Table.read("meteors_database.tab", format="ascii")
meteor_mask = load_meteor_pos["Fecha"] == date
plt.plot(load_meteor_pos["Longitud"][meteor_mask], load_meteor_pos["Latitud"][meteor_mask], "mo")


# Develop dictionaries with events positions, their respective stations.

dates_stations_dict = {"2019-05-23":("KVTX", "TNCU", "UAGU"), "2019-07-18":("KVTX", "MDO1", "TNCU", "UAGU"), "2019-08-10":("KVTX", "TNCU", "UAGU", "UCOE", "UGEO"), "2019-10-03":("KVTX", "UAGU", "UXAL"), "2019-10-09":("GUAX", "TNBA", "TNHM","TNMS"), "2019-11-17":("GUAX", "TNBA", "PTEX", "TSFX", "TNPP", "P001", "BAR1", "WWMT"), "2019-11-16":("KVTX", "P807", "SG33", "TNCU", "UAGU"), "2019-11-19": ("CN23", "GCFS", "UNPM"), "2019-11-26":("TNBA", "TNHM", "TNMS", "UAGU", "YESX"), "2019-12-04":("GMPK", "IAGX", "P014", "PLPX", "TNHM", "TNPP"), "2019-12-15":("GUAX", "TNBA", "TNCU", "TNHM", "TNMS", "UAGU", "YESX"), "2019-12-29":("BAR1", "BLYT", "GUAX", "P014", "PTEX", "TNBA", "TNHM", "USMX", "WWMT"), "2020-01-03":("BAR1", "BLYT", "GUAX", "P014", "PTEX", "TNBA", "TNHM", "USMX", "WWMT"), "2020-01-06":("P014", "P807", "RG07", "TNCU", "USMX", "YESX"), "2020-01-15":("CN23", "CN25", "OXTH", "TNAT", "TNNX", "UNPM", "UXAL"), "2020-02-12":("CN23", "GUAT", "UNPM", "OXTH", "TNNX", "UXAL"), "2020-03-03":("TNAM", "TNCC", "TNCM", "TNMS"), "2020-03-31":("GUAX", "TNCU", "TNHM", "TNPP"), "2020-04-08":("KVTX", "MGW3", "UHSL", "UNPM", "UXAL"), "2020-04-18":("KVTX", "SG33", "TNCU", "TNHM", "UAGU", "YESX"), "2020-04-20":("KVTX", "MGO5", "MGW3", "P807", "WEPD", "WMOK"), "2020-04-25":("P001", "P014", "RG06", "TNPP", "USMX"), "2020-04-28":("TNCM", "TNNP", "YESX"), "2020-05-08":("KVTX", "MGW3", "UNPM", "UXAL"), "2020-07-15":("GUAX", "INEG", "TNAM", "TNCU", "TNHM", "TNMS", "YESX"), "2020-08-07":("INEG", "KVTX", "MGO5", "SG33", "TNCU", "USMX"), "2020-09-13":("GUAX", "PTEX", "TNHM", "TSFX"), "2020-09-30":("GUAX", "INEG", "TNAM", "TNCU", "TNHM", "USMX"), "2020-11-16":("INEG", "TNAM", "TNCN", "TNGF", "UCOE"), "2020-11-17":("INEG", "P807", "TNAM", "TNCU", "UCOE"), "2020-12-19":("INEG", "TNAM", "TNCU", "UCOE", "UHWL", "UXAL"), "2020-12-23":("GUAX", "PTEX", "TNHM", "TSFX"), "2020-12-29":("TNCN", "TNGF", "TNNX", "TNSJ"), "2021-03-31":("OXUM", "TGMX", "TNNX", "UXAL")}

# Plot stations positions

stations = dates_stations_dict[date]

for station in stations:
    mask = stations_pos["Site"] == station
    plt.plot(stations_pos["Longitude"][mask], stations_pos["Latitude"][mask], "go")
    plt.annotate(station, (stations_pos["Longitude"][mask][0], stations_pos["Latitude"][mask][0]),
			  textcoords="offset points", color="w", xytext=(10, 10), ha="center", bbox=dict(boxstyle="round", pad=0.5, fc="g", alpha=0.7))


# Plot bolide trajectory

GLM16_file = open(directory+"GLM-16-data.csv")
GLM17_file = open(directory+"GLM-17-data.csv")

for i in range(10): # skip unneeded data
    GLM16_file.readline()
    GLM17_file.readline()

GLM16_data = GLM16_file.readlines()
GLM17_data = GLM17_file.readlines()
GLM16_table = Table.read(GLM16_data, format="ascii")
GLM17_table = Table.read(GLM17_data, format="ascii")

f1_longitude, f1_latitude = GLM16_table["longitude"], GLM16_table["latitude"]
f2_longitude, f2_latitude = GLM17_table["longitude"], GLM17_table["latitude"]


fit_coord1 = np.polyfit(f1_longitude, f1_latitude, 1)
fit_coord2 = np.polyfit(f2_longitude, f2_latitude, 1)


poly1 = np.poly1d(fit_coord1)
poly2 = np.poly1d(fit_coord2)

if f1_longitude[-1] > f1_longitude[0]:
    step1 = -2
else:
    step1 = 2

if f2_longitude[-1] > f2_longitude[0]:
    step2 = -2
else:
    step2 = 2
#step = 0.5*(step1+step2)

plt.plot([f1_longitude[0]+step1, f1_longitude[-1]], [poly1(f1_longitude[0]+step1), poly1(f1_longitude[-1])], "r--")
plt.plot([f2_longitude[0]+step2, f2_longitude[-1]], [poly2(f2_longitude[0]+step2), poly2(f2_longitude[-1])], "r--")

#plt.plot([0.5*(f1_longitude[0]+f2_longitude[0])+step, 0.5*(f1_longitude[-1]+f2_longitude[-1])], [0.5*(poly1(f1_longitude[0]-2)+poly2(f2_longitude[0]+step)), 0.5*(poly1(f1_longitude[-1])+poly2(f2_longitude[-1]))], "k", lw=2)
plt.annotate(date, (0.5*(f1_longitude[-1]+f2_longitude[-1]), 0.5*(poly1(f1_longitude[-1])+poly2(f2_longitude[-1]))),
			  textcoords="offset points", color="w", xytext=(10, 10), ha="center", bbox=dict(boxstyle="round", pad=0.5, fc="r", alpha=0.7))

ax = plt.gca()
ax.set_aspect("equal", adjustable="box")
outfolder = "./stations_maps/"
plt.savefig(outfolder+date+"-stations.pdf")
