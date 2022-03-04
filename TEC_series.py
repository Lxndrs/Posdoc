import pandas as pd
import argparse 
import numpy as np
import matplotlib.pyplot as plt
import glob
from statistics import mode
import seaborn as sns
parser = argparse.ArgumentParser(
	   description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
				      help='Choose date. Format: yyyy-mm-dd')


#parser.add_argument("--filterhigh", type=float, default=30, 
#			    help="Filter the upper time interval of the TEC curve (in minutes)")
#parser.add_argument("--filterlow", type=float, default=10, 
#			    help="Filter the lower time interval of the TEC curve (in minutes)")
parser.add_argument("--previous", action="store_true", help="Choose previous day data")
parser.add_argument("--sets", type=str, help="Choose between sets of stations")
cmd_args = parser.parse_args()
date = cmd_args.date
#filter_high = cmd_args.filterhigh
#filter_low = cmd_args.filterlow
previous = cmd_args.previous
stations_set = cmd_args.sets
set_folder = "/set"+stations_set
directory = "./data/"+date+set_folder
outdir = "./TEC_series/"+date+set_folder
infile = glob.glob(directory+"/*.csv")[0]
Date, time, residual = infile.split("_")
if previous:
    infile = glob.glob(directory+"/previous/*.csv")[0]
    outdir = "./TEC_series/"+date+set_folder+"/previous"
# ****************************** Read data file ************************************************

csv_file = pd.read_csv(infile)
station_array = np.unique(csv_file["Station"])
time_hours = float(time)
fh = 130
fl = -60
ws = fh-fl
sns.set_style("whitegrid")
#for ws in window_size:
#    win_mask = csv_file["WSize"] == ws
for s in station_array:
    s_mask = csv_file["Station"] == s
    prn_array = np.unique(csv_file["PRN"][s_mask])
    for p in prn_array:
        prn_mask = csv_file["PRN"] == p
        frequency = round(mode(np.gradient(csv_file["Time"][s_mask & prn_mask]*3600.)))
        print("Station:{}, PRN:{}, Frequency:{}".format(s, p, frequency))
        plt.plot(csv_file["Time"][prn_mask & s_mask], csv_file["vTEC"][prn_mask & s_mask])
        plt.axvline(time_hours, ls="--", c="k")
        plt.axhline(0, ls="--", c="k")
        plt.ylim(-0.3, 0.5)
        plt.xlabel("UT (hours)")
        plt.ylabel("vTEC (TECU)")
        plt.title(s.upper()+"-PRN{}".format(p))
        plt.savefig(outdir+"/{}-TEC_series_PRN_{}.pdf".format(s, p))
        plt.clf()
