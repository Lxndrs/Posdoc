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

parser.add_argument("--formato", type=str, default="pdf", choices=("pdf", "png", "jpg"), 
                                help="Choose output format")
parser.add_argument("--starttime", type=str, default="00:00:00", help="Select the event start time")
parser.add_argument("--filterhigh", type=float, default=30, 
                     help="Filter the upper time interval of the TEC curve (in ninutes)")
parser.add_argument("--filterlow", type=float, default=10, 
                     help="Filter the lower time interval of the TEC curve (in ninutes)")

cmd_args = parser.parse_args()
date = cmd_args.date
formato = cmd_args.formato
start_time = cmd_args.starttime
filter_high = cmd_args.filterhigh
filter_low = cmd_args.filterlow

# ****************************** Read data file ************************************************

f = open("./data/2019-02-01/unpm032-2019-02-01.Cmn") # Test file until the program works
for i in range(4):
    f.readline()

raw_data = f.readlines()

# ****************************** Extract relevant information ***********************************

data = Table.read(raw_data, format="ascii")
hr, minute, sec = start_time.split(":")
time_hours = float(hr)+float(minute)/60. + float(sec)/3600.
time = data["Time"]
vTEC = data["Vtec"]
PRN = data["PRN"]
time_corrector = time < 0
time[time_corrector] = time[time_corrector] + 24.0
tau_0, zeta_0 = 2.0, 40.0
prn_array = np.unique(PRN)
for p in prn_array:
    PRN_mask = PRN == p
    selected_time = time[PRN_mask]
    selected_TEC = vTEC[PRN_mask]
    time_mask = (selected_time < time_hours + filter_high/60.) & (selected_time > time_hours-filter_low/60.)

# ******************************* Split the signal into fragments *********************************

    S_time, S_TEC = split_PRN(selected_time, selected_TEC)
#test plot
    for s, tec in zip(S_time, S_TEC):
        plt.plot(s, tec, "b.")
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
#        fdp_vt = interp1d(brc_t, brc_vt) # The free depletion signal may be an interpolation of the BRC
        fdp_x, fdp_y = free_dep_signal(brc_t, brc_vt, s, tec)  
        plt.plot(brc_t, brc_vt, "r*")
        plt.plot(fdp_x, fdp_y)
#        plt.plot(s, fdp_vt(s))
# ********************************** Get trend using Savitzky-Golay filter ***********************
#        residuals = 1.0
#        ord=1
#        while(residuals > 0.15)|(ord <=10): # Use SG filter for the resulting fragments. Increase order of
                                        # fit until residuals become small enough
#            y_trend = savitsky(tec, window_length=get_window_size(len(tec)), polyorder=ord)
#            residuals = np.sum((y_trend-tec)**2)/len(tec)
#            ord=ord+1
#        plt.plot(s, y_trend, "r--")
#        print(residuals)
    plt.savefig("./TEC_tests/test_detrend_{}.pdf".format(p))
    plt.clf()
#plt.savefig("./TEC_tests/test_detrend_test.pdf")

 

#    print(p, len(X))
#if len(X) < 20: # TEC series with too few data will be considered "empty"
#    continue







#, vt = X*tau_0, Y*zeta_0




 
# ********************************** Plot test graph *********************************************
#    print(p, np.corrcoef(t, vt)[0][1])
#plt.plot(t, vt, "g.")


#plt.plot(t, y_trend, "m--")
#    plt.plot(X, Y, "g.")
#    plt.plot(BRC_x, BRC_y, "r*")
#    plt.fill_between(np.array(BRC_x), np.array(BRC_y)+1, np.array(BRC_y)-3, alpha=0.1) 
#    plt.axvline(time_hours, ls="--", c="k")
