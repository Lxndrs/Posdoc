import numpy as np
from astropy.table import Table
import argparse 
import matplotlib.pyplot as plt
import seaborn as sns

# Script goal: Plot Dst index vs time
#                                      WDC for Geomagnetism, Kyoto
#                                Hourly Equatorial Dst Values (REAL-TIME)  

# Read arguments from Terminal and open file

parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
			       help='Choose date. Format: yyyy-mm-dd')

cmd_args = parser.parse_args()
date = cmd_args.date
year, month, day = date.split("-")
datafile = "Dst-index-"+month+"-"+year+".txt"
directory = "./Dst index/" 
f = open(directory+datafile, "r")

# Extract information from file

h1 = f.readline()
print(h1)
h2 = f.readline()
print(h2)
h3 = f.readline()
print(h3)
f.readline()
time_str = f.readline()
f.readline()
data = f.readlines()

Tabdata= Table.read(data, format="ascii")
Dst_index = []

mask = Tabdata["col1"]==int(day)
maskp = Tabdata["col1"]==int(day) - 1
# Set the time array for plotting

time = time_str.split()
for i in range(len(time)):
    time[i] = int(time[i])
#time=timep+24
#time_c = np.array([timep, time])
#time_c = np.reshape()
# Set the Dst index array for plotting
for k in Tabdata.keys()[1:]:
    Dst_index.append(Tabdata[k][mask][0])

# Plot Dst-index vs time

sns.set_style("whitegrid")
plt.plot(time, Dst_index)
plt.xlim(0, 24)
plt.ylim(-500, 100)
plt.title("DST (Final)")
plt.text(0, 85, day+" "+h3.lstrip())
plt.text(15, 85, h1.lstrip())
plt.xlabel("Time (UT, hours)")
plt.ylabel("Dst Index (nT)")
outfile = datafile.replace(".txt", ".pdf")
plt.savefig(directory+outfile)
