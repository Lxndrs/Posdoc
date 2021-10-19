import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.table import Table
import argparse
import glob

# Goal of this program: Obtain TEC vs time plots and analyze them in order to check if we can derive some kind
# of parameter to quantify the TEC perturbations

# We must recycle a lot of the code lines from plot_vTEC.py

parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
				   help='Choose date. Format: yyyy-mm-dd')




cmd_args = parser.parse_args()
date = cmd_args.date


directory = "./data/"+date
p_directory = directory + "/previous/"
n_directory = directory+ "/next/"

# set figure style
sns.set_style("whitegrid") 

# Load RINEX capabilities


std_files = glob.glob(directory+"/*.Std")
load_std = [Table.read(std_files[i], format="ascii") for i in range(len(std_files))]


std_p = glob.glob(p_directory+"*.Std")
std_n = glob.glob(n_directory+"*.Std")

load_std_p = [Table.read(std_p[i], format="ascii") for i in range(len(std_p))]
load_std_n = [Table.read(std_n[i], format="ascii") for i in range(len(std_n))]

# Load backgrond median



# Get the data and plot

fig = plt.figure()
j=1
for f, fp, fn in zip(load_std, load_std_p, load_std_n):
    std_time = f["col1"]
    std_time_p = fp["col1"]
    std_time_n = fn["col1"]
    std_TEC = f["col2"]
    std_TEC_p = fp["col2"]
    std_TEC_n = fn["col2"]

    for i in range(len(std_TEC)):
        if std_TEC[i] == "-":
            std_TEC[i]=np.nan
    for i in range(len(std_TEC_p)):
        if std_TEC_p[i] == "-":
            std_TEC_p[i]=np.nan
    for i in range(len(std_TEC_n)):
        if std_TEC_n[i] == "-":
            std_TEC_n[i]=np.nan

    avg_TEC = np.array([std_TEC_p, std_TEC, std_TEC_n])
    avg_TEC = np.reshape(avg_TEC, len(std_TEC_p)+len(std_TEC)+len(std_TEC_n))

    avg_time = np.array([std_time_p, std_time+24.0, std_time_n+48.0])
    avg_time = np.reshape(avg_time, len(std_time_p)+len(std_time)+len(std_time_n))
    plt.plot(avg_time, avg_TEC)   
    ax.set_ylim(0, 40)

fig.savefig("TEC.pdf")
