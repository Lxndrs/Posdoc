# This program is designed to estimate the standard deviation of
# the Rate of TEC (ROTI) from 1 Hz RINEX data for every PRN in a
# *.Cmn file (output of the GPS GOPI software where we get TEC 
# calculations from Hatanaka files)

import argparse
import numpy as np
from astropy.table import Table
import glob

# Import parameters from command line -----------------------------------------

parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
			       help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--station", type=str, default="kvtx", 
                                help="Choose station to work")



cmd_args = parser.parse_args()
date = cmd_args.date
station = cmd_args.station

# Read file to work -----------------------------------------------------------

directory = "./data/"+date+"/high-rate/"
pathfile = glob.glob(directory+station+"*.Cmn")
f = open(pathfile[0], "r") #pathfile is an array with one element in this case
f.readline()
f.readline()
f.readline()
f.readline() # The first 4 lines do not contain useful information
data = f.readlines()
Tabdata = Table.read(data, format="ascii")
time = Tabdata["Time"]
time_mask = time < 0
time[time_mask] = time[time_mask] + 24.0 # Fix time values below 0

# Get an array with the available PRNs, set output arrays and loop ----------------

PRNs = np.unique(Tabdata["PRN"])
N = 60 #The ROTI will be computed on a 60 seconds time interval
ROTI = []
outfile = "./data/"+date+"/ROTI-compute_"+station+"_"+date+".tab"
for p in PRNs:
    PRN_mask = Tabdata["PRN"] == p # select only the data of the same PRN
    f_time = np.round(time[PRN_mask]*3600) # time will be measured 
                                           # in seconds and avoid decimals
    f_stec = Tabdata["Stec"][PRN_mask] # the ROTI will be computed using sTEC data
#    t_min, t_max = min(f_time), min(f_time)+N #Define time window where ROTI is 
                                              #computed
    ROT = np.gradient(f_stec)/np.gradient(f_time) #Compute ROT for the whole PRN
    for i in range(N):
        ROTI.append(np.nan) # Before the first minute we can't compute ROTI
    for i in range(len(f_time)-N):
        r_mask = (f_time <= f_time[i+N]) & (f_time >= f_time[i])
        roti = np.std(ROT[r_mask])
        ROTI.append(roti) 

#    while(t_max <= f_time[-1]):#when the upper limit reaches the end, stops this loop
#        r_mask = (f_time <= t_max) & (f_time >= t_min) # This loop won't work because the 1 sec
#        roti = np.std(ROT[r_mask])                     # time interval doesn't follows all the time
#        ROTI.append(roti)
#        t_min = t_min+1
#        t_max = t_max+1 # Shifts the interval time window in 1 second

# Save to a file -------------------------------------------------------------------


outtab = Table([time, Tabdata["PRN"], Tabdata["Az"], Tabdata["Ele"], Tabdata["Lat"], Tabdata["Lon"], Tabdata["Stec"], Tabdata["Vtec"], ROTI], names=("Time", "PRN", "Az", "Ele", "Lat", "Lon", "Stec", "Vtec", "ROTI"))
outtab.write(outfile, format="ascii")
