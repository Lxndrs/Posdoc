import numpy as np
from astropy.table import Table
import glob

# The goal of this program is to estimate the total energy and peak energy
# from GLM data

# Methodology: sum all the values from the energy column of GLM data and also
# get the maximum value. Both should be outputs of the program

# Get a list of the paths of al GLM data

dates = ["2019-02-01", "2019-05-23", "2019-07-18", "2019-08-10", "2019-10-03", "2019-10-09", "2019-11-16", "2019-11-17", "2019-11-19", "2019-11-26", "2019-12-04", "2019-12-15", "2019-12-29", "2020-01-03", "2020-01-06", "2020-01-15", "2020-02-12", "2020-03-03", "2020-03-31", "2020-04-08", "2020-04-18", "2020-04-20", "2020-04-25", "2020-04-28", "2020-05-08", "2020-07-15", "2020-08-07", "2020-09-13","2020-09-30", "2020-11-16", "2020-11-17", "2020-12-19","2020-12-23", "2020-12-29", "2021-03-31"]

G16_dirs = [] 
G17_dirs = []
root_folder = "./data/"
common_folder = "/GLM/"
for date in dates:
    g16 = glob.glob(root_folder+date+common_folder+"*16*")
    g17 = glob.glob(root_folder+date+common_folder+"*17*")
    G16_dirs.append(g16[0])
    G17_dirs.append(g17[0])
# initialize output arrays
g16_total_energy_array = []
g17_total_energy_array = []


#initialize loop for data adquisition

for g16, g17 in zip(G16_dirs, G17_dirs):
    data16 = open(g16, "r")
    data17 = open(g17, "r")

    for i in range(10): # skip unneded data
        data16.readline()
        data17.readline()

    # gather the table with meteor info
    g16data = data16.readlines()
    g17data = data17.readlines()

    g16table = Table.read(g16data, format="ascii")
    g17table = Table.read(g17data, format="ascii")

    # estimate peak energy and total energy for each satellite
    g16_total_energy = np.sum(g16table["energy (joules)"])
    g16_total_energy_array.append(g16_total_energy)
    g17_total_energy = np.sum(g17table["energy (joules)"])
    g17_total_energy_array.append(g17_total_energy)


print("Date (yyyy-mm-dd)", "GLM-16 total energy (joules)","GLM-17 total energy (joules)")
for date, g16_total, g17_total in zip(dates, g16_total_energy_array, g17_total_energy_array):
    print(date, g16_total, g17_total)
