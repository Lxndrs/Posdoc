# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# By Alejandro Tarango-Yong
# This script is the core of the TEC maps creator. The idea is to load the already created
# files for the compilations of the detrended files and collect the data from the specific instant
# chosen in the command line.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Import python libraries

import numpy as np
import pandas as pd
import argparse
from astropy.table import Table

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extract the data from the command line

parser = argparse.ArgumentParser(description = """Extract command line arguments""")
parser.add_argument("--time", type=str, default="21:25:48", help="Instant where map is taken")
parser.add_argument("--TEC", type=str, default="vTEC", choices=("vTEC", "sTEC"), 
                    help="Choose type of TEC data")

cmd_args = parser.parse_args()
time_snap = cmd_args.time
TEC_choice = cmd_args.TEC
TEC_dict = {"vTEC":"TEC", "sTEC":"sTEC"}
hrs, mins, secs = time_snap.split(":")
time_snap_dec = float(hrs) + float(mins)/60. + float(secs)/3600.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Read csv files for the 5 sets of stations (We work with what we have)
#
base_folder = "./data/2019-06-22/"
base_file = "2019-06-22_21.43_detrended"
set1 = pd.read_csv("{}set1/{}-{}.csv".format(base_folder, base_file, TEC_dict[TEC_choice]))
set2 = pd.read_csv("{}set2/{}-{}.csv".format(base_folder, base_file, TEC_dict[TEC_choice]))
set3 = pd.read_csv("{}set3/{}-{}.csv".format(base_folder, base_file, TEC_dict[TEC_choice]))
set4 = pd.read_csv("{}set4/{}-{}.csv".format(base_folder, base_file, TEC_dict[TEC_choice]))
set5 = pd.read_csv("{}set5/{}-{}.csv".format(base_folder, base_file, TEC_dict[TEC_choice]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extract and merge relevant data
#
extract_if = (set1["Time"] < time_snap_dec + 7.5/3600.) & (set1["Time"] > time_snap_dec - 7.5/3600)
# Only extract the data within an interval of 15 seconds where time_snap is the center in case
# that the chosen time did not match with the GPS grid.

lat1 = np.array(set1["Latitude"][extract_if])
lat2 = np.array(set2["Latitude"][extract_if])
lat3 = np.array(set3["Latitude"][extract_if])
lat4 = np.array(set4["Latitude"][extract_if])
lat5 = np.array(set5["Latitude"][extract_if])
latitude = np.concatenate((lat1, lat2, lat3, lat4, lat5))

lon1 = np.array(set1["Longitude"][extract_if])
lon2 = np.array(set2["Longitude"][extract_if])
lon3 = np.array(set3["Longitude"][extract_if])
lon4 = np.array(set4["Longitude"][extract_if])
lon5 = np.array(set5["Longitude"][extract_if])
longitude = np.concatenate((lon1, lon2, lon3, lon4, lon5))

tec1 = np.array(set1[TEC_choice][extract_if])
tec2 = np.array(set2[TEC_choice][extract_if])
tec3 = np.array(set3[TEC_choice][extract_if])
tec4 = np.array(set4[TEC_choice][extract_if])
tec5 = np.array(set5[TEC_choice][extract_if])
TEC = np.concatenate((tec1, tec2, tec3, tec4, tec5))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Save everything in a single file

output_file = "{}{}_snapshot_{}.tab".format(base_folder, TEC_choice, time_snap)
output_table = Table([latitude, longitude, TEC], names = ("Latitude", "Longitude", TEC_choice))
output_table.write(output_file, format="ascii", overwrite=True)
print("{} map file completed, time = {}".format(TEC_choice, time_snap))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

