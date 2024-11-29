# This program is designed to estimate the standard deviation of
# the Rate of TEC (ROTI) from low resolution RINEX data for every PRN in a
# *.Cmn file (output of the GPS GOPI software where we get TEC 
# calculations from Hatanaka files)

import argparse
import numpy as np
from astropy.table import Table
import glob
import pandas as pd

# Import parameters from command line -----------------------------------------

parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
			       help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--station", type=str, default="kvtx", 
                                help="Choose station to work")

parser.add_argument("--sets", type=str, default="1", help="Choose data set")
parser.add_argument("--station", type=str, default="bara", help="Choose station")

cmd_args = parser.parse_args()
date = cmd_args.date
station = cmd_args.station
set_folder= cmd_args.sets
station = cmd_args.station

# Read file to work -----------------------------------------------------------

directory = "./data/{}/set{}/".format(date, set_folder)
pathfile = glob.glob(directory+"{}*-{}.Cmn".format(date, station))
f = open(pathfile[0], "r")
f.readline()
f.readline()
f.readline()
f.readline() # The first 4 lines do not contain useful information
data = f.readlines()
Tabdata = Table.read(data, format="ascii")
time = Tabdata["Time"]
tabPRN = Tabdata["PRN"]
tabsTEC = Tabdata["Stec"]
tabLat = Tabdata["Lat"]
tabLon = Tabdata["Lon"]
outfile = "./data/"+date+"/ROTI-compute_"+station+"_"+date+".csv"
time_mask = time < 0
time[time_mask] = time[time_mask] + 24.0 # Fix time values below 0

# Get an array with the available PRNs, set output arrays and loop ----------------

PRNs = np.unique(tabPRN)
time_window = 60 #The ROTI will be computed on a 60 seconds time interval
ROTI = np.nan*np.zeros_like(time)

for p in PRNs:
    PRN_mask = tabPRN == p # select only the data of the same PRN
    f_time = np.round(time[PRN_mask]*3600) # time will be measured 
                                           # in seconds and avoid decimals
    f_stec = tabsTEC[PRN_mask] # the ROTI will be computed using sTEC data
#    t_min, t_max = min(f_time), min(f_time)+N #Define time window where ROTI is 
                                              #computed
    ROT = np.gradient(f_stec)/np.gradient(f_time) #Compute ROT for the whole PRN

#    for i in range(N):
#        ROTI.append(np.nan) # Before the first minute we can't compute ROTI
#    for i in range(len(f_time)-N):
#        r_mask = (f_time <= f_time[i+N]) & (f_time >= f_time[i])
#        roti = np.std(ROT[r_mask])
#        ROTI.append(roti) 


# Save to a file -------------------------------------------------------------------


outtab = Table([time, tabPRN, tabLat, tabLon, 
         tabsTEC, ROTI], names=("Time", "PRN", "Lat", "Lon", "Stec", "ROTI"))
outtab.write(outfile, format="csv", overwrite=True)
