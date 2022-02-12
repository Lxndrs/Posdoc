import numpy as np
from astropy.table import Table
import argparse
import glob 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter as savitsky
from scipy.stats import mode

# The goal of this progrma is to get detrended time series TEC curves using the method
# Described by Pradipta et al. 2015. We expect that the resulting detrended TEC data show
# ionospheric perturbations in our dataset.


def next_contact_point(x, y, x0, y0):
    """
    Estimate the next contact point into the Barrell roll curve (BRC)
    inputs:
    x: scaled time.
    y: scaled vTEC
    x0: x pivot point
    y0: y pivot point
    output:
    (xf, yf): coordinates of the next contact point
    """
    # Define region of interest
    R0 = 1 # Radius of the barrell roll. Unitary radius works pretty well
    ROI = (x > x0) & (x < x0 + 2*R0)
    # Delta_x and delta_y are the separation between the elements which belong to the ROI and (x0, y0)
    delta_x = x[ROI]-x0
    delta_y = y[ROI]-y0

    try:
	  # calculating important angles
        theta = np.arctan2(delta_y, delta_x)
        cos_alpha = np.sqrt(delta_x**2 + delta_y**2)/(2*R0) 
        delta = np.arcsin(cos_alpha) - theta
	  # Selecting the Next Contact Point (NCP)
        NCP = delta == min(delta) # The next contact point has the smallest angular distance delta
        xf, yf = x0+delta_x[NCP][0], y0+delta_y[NCP][0]
    except ValueError: # this happens because ROI is empty
        xf, yf = x[~ROI & (x>x0)][0], y[~ROI & (x>x0)][0] #The element we will use as next contact point
							    # is the first outside the ROI and greater than x0
    return xf, yf

def free_dep_signal(brc_x, brc_y, x, y):
    """
    Create the free depletion signal.
    Inputs:
    (brc_x, brc_y): Barrel roll curve coordinates
    (x, y): TEC curve 
    output:
    (xf, yf): free depletion signal
    """
    delta_1, delta_2 = 1, 3
    BRC_plus = brc_y + delta_1
    BRC_minus = brc_y - delta_2
    int_bplus = interp1d(brc_x, BRC_plus)
    int_bminus = interp1d(brc_x, BRC_minus)
    mask = (y < int_bplus(x)) & (y> int_bminus(x))
    xf, yf = x[mask], y[mask]
    return xf, yf

# ************************** Get window size for Savitzky-Golay Filter ************************

def get_window_size(x):
    """
    Get window size for the Savitzky-Golay filter. Must be an *odd* integer.
    This module basically substracts 1 if x is an even number and keeps the number if x is odd
    the output number is integer
    """
    return int(2*np.ceil(x/2.)-1)

def discontinuity(x):
    """
    Find discontinuities in any array
    """
    x_gradient = np.gradient(x)
    dx_mode = mode(x_gradient)[0][0]
    for dx in x_gradient:
        if dx > 2*dx_mode:
            flag = True
            break
        else:
            flag = False
    return flag

# **************************** Split discontine PRN *******************************************

def split_PRN(t, TEC):
    """
    Split a discontinuos PRN into pieces and return the biggest ones
    inputs:
    t, TEC: arrays to be split. t --> time , TEC --> vTEC
    outputs:
    output_t, output_TEC --> 2D arrays which contain the fragmented curve
    """

    index_discontinuity =[]
    gradient_t = np.gradient(t)
    for i, dt in enumerate(gradient_t):
        if dt > 0.01: # 0.01 is like 10 times the regular GPS frequency
            index_discontinuity.append(i) # collect the indices where time gradient is big enough
    split_t = np.split(t, index_discontinuity)
    split_TEC = np.split(TEC, index_discontinuity)
    output_t = []
    output_TEC =[]
    for s, tec in zip(split_t, split_TEC):
        if len(s) > 20: #if the subarray contain too few elements will be discarded
            output_t.append(s)
            output_TEC.append(tec)
    return output_t, output_TEC

# *************************** Read inputs from command line ***********************************

parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
				 help='Choose date. Format: yyyy-mm-dd')
parser.add_argument("--starttime", type=str, default="00:00:00", help="Select the event start time")
parser.add_argument("--previous", action="store_true", help="Choose previous day data")
parser.add_argument("--sets", type=str, help="Choose between sets of stations")

cmd_args = parser.parse_args()
date = cmd_args.date
start_time = cmd_args.starttime
previous = cmd_args.previous
filter_low = np.array([-15, -35, -60])
filter_high = np.array([40, 80, 130])
stations_set = cmd_args.sets
set_folder = "/set"+stations_set
# ****************************** Read data file ************************************************

directory = "./data/"+date+set_folder
if previous:
    directory = "./data/"+date+set_folder+"/previous"
files = glob.glob(directory+"/*.Cmn")

# ******************************* Set final output lists ******************************************

stations = []
final_time = []
final_PRN = []
final_lat = []
final_lon = []
final_vTEC = []
window_size = []

# Use the longest window size to keep most frequencies inside COI
#for fl, fh in zip(filter_low, filter_high): # Now the window size is variable
fh = 130
fl = -60
ws = fh-fl
for File in files:
    f = open(File, "r")
    for i in range(4):
        f.readline()

    raw_data = f.readlines()

# ****************************** Extract relevant information ***********************************
    station = File[-22:-18] # Extract station name from file name(This should work in both
			    # windows and linux)
    data = Table.read(raw_data, format="ascii")
    hr, minute, sec = start_time.split(":")
    time_hours = float(hr)+float(minute)/60. + float(sec)/3600.
    time = data["Time"]
    vTEC = data["Vtec"]
    PRN = data["PRN"]
    latitude = data["Lat"]
    longitude = data["Lon"] -360.
    time_corrector = time < 0
    time[time_corrector] = time[time_corrector] + 24.0
    tau_0, zeta_0 = 2.0, 40.0
    prn_array = np.unique(PRN)
    time_mask = (time < time_hours + fh/60.) & (time > time_hours-fl/60.)

# ******************************* Start loop over PRNs ********************************************

    for p in prn_array:
        PRN_mask = PRN == p
        filtered_time = time[PRN_mask & time_mask]
        filtered_TEC = vTEC[PRN_mask & time_mask]
        filtered_lat = latitude[PRN_mask & time_mask]
        filtered_lon = longitude[PRN_mask & time_mask]
        if len(filtered_time) < 20: # TEC series with too few data will be considered "empty"
            continue

# ******************************* Split the signal into fragments *********************************

        S_time, S_TEC = split_PRN(filtered_time, filtered_TEC)
        S_time, S_lat = split_PRN(filtered_time, filtered_lat)
        S_time, S_lon = split_PRN(filtered_time, filtered_lon)
        for s, tec, la, lo in zip(S_time, S_TEC, S_lat, S_lon):
# ******************************* Start making BRC ***********************************************
            X, Y = s/tau_0, tec/zeta_0
            x_0, y_0 = X[0], Y[0]
            BRC_x = [x_0]
            BRC_y = [y_0]
            while(x_0 < X[-1]):
                xn, yn = next_contact_point(X, Y, x_0, y_0)
                BRC_x.append(xn)
                BRC_y.append(yn)
                x_0, y_0 = xn, yn

# ********************************** Making free deptetion signal ********************************
# Return to the TEC-time space
            brc_t, brc_vt = np.array(BRC_x)*tau_0, np.array(BRC_y)*zeta_0
            fdp_x, fdp_y = free_dep_signal(brc_t, brc_vt, s, tec)  
# ********************************** Get trend using Savitzky-Golay filter **************************
            y_trend = savitsky(tec, window_length=get_window_size(len(tec)), polyorder=7)
            residuals = np.sum((y_trend-tec)**2)/len(tec)
            plt.plot(s, y_trend, "r--", label="Trend")
# *********************************** Detrend signal and plot undetrended signal and trend************
            det_signal = tec - y_trend
            plt.plot(s, tec, label="original signal")
            plt.axvline(time_hours, c="k", ls="--")
            plt.xlabel("UT (hours)")
            plt.ylabel("vTEC (TECU)")
            plt.legend()
            plt.title("{} PRN{}. Date: {}".format(station, p, date))
            plt.savefig(directory+"/trends/{}-{}-{}.pdf".format(station, p, ws))
            plt.clf()
            for ft, flat, flon, d in zip(s, la, lo, det_signal): # Fill output lists
                final_vTEC.append(d)
                final_time.append(ft)          
                final_PRN.append(p)
                final_lat.append(flat) 
                final_lon.append(flon) 
                stations.append(station)
    #            window_size.append(ws)
# ***************************** Mask positions and save data to file *******************************


output_table = Table([stations, final_time, final_PRN, final_lat, final_lon, final_vTEC], names=("Station", "Time", "PRN", "Latitude", "Longitude", "vTEC"))
outfile = directory+"/"+date+"_"+str(time_hours)+"_detrended-TEC.csv" #File.replace("Cmn", "csv")
output_table.write(outfile, format="csv", overwrite=True)
