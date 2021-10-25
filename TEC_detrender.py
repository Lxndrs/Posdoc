import numpy as np
from astropy.table import Table
import argparse
import glob 
import matplotlib.pyplot as plt

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
    # calculating important angles
    theta = np.arctan2(delta_y, delta_x)
    cos_alpha = np.sqrt(delta_x**2 + delta_y**2)/(2*R0) 
    delta = np.arcsin(cos_alpha) - theta
    # Selecting the Next Contact Point (NCP)
    try:
        NCP = delta == min(delta) # The next contact point has the smallest angular distance delta
        return x0+delta_x[NCP][0], y0+delta_y[NCP][0]
    except ValueError: # IF the ROI is empty use the first point outside the ROI (besides (x0, y0))
        return x[~ROI][1], y[~ROI][1]

# *************************** Read inputs from command line ***********************************

parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
				 help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--formato", type=str, default="pdf", choices=("pdf", "png", "jpg"), 
				  help="Choose output format")

cmd_args = parser.parse_args()
date = cmd_args.date
formato = cmd_args.formato

# ****************************** Read data file ************************************************


f = open("./data/2019-02-01/unpm032-2019-02-01.Cmn") # Test file until the program works

for i in range(4):
    f.readline()

raw_data = f.readlines()

# ****************************** Extract relevant information ***********************************

data = Table.read(raw_data, format="ascii")
time = data["Time"]
vTEC = data["Vtec"]
PRN = data["PRN"]
time_corrector = time < 0
time[time_corrector] = time[time_corrector] + 24.0
PRN_mask = PRN ==1
# ******************************* Start making BRC **********************************************
tau_0, zeta_0 = 2.0, 40.0
X, Y = time[PRN_mask]/tau_0, vTEC[PRN_mask]/zeta_0
x_0, y_0 = X[0], Y[0]
BRC_x = []
BRC_y = []
j=0
while((x_0 < X[-1]) | j==500):
    xn, yn = next_contact_point(X, Y, x_0, y_0)
    BRC_x.append(xn)
    BRC_y.append(yn)
    x_0, y_0 = xn, yn
    j=j+1 #stop the while cycle in case of it fails
# ********************************** Plot test graph *********************************************

plt.plot(time[PRN_mask], vTEC[PRN_mask], "g.")
plt.plot(np.array(BRC_x)*tau_0, np.array(BRC_y)*zeta_0, "r*")
plt.savefig("test_detrend.pdf")
