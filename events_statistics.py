import numpy as np
from astropy.table import Table
import statistics as stats
import matplotlib.pyplot as plt

# Code goal: Extract basic statistics of events sample

# Step 1: Read data from table

tab_data = Table.read("meteors_database.tab", format="ascii")

duration = tab_data["dt"]
s_dur = tab_data["s_dt"]
lat = tab_data["Latitud"]
lon = tab_data["Longitud"]
s_lat = tab_data["sig_latitude"]
s_lon = tab_data["sig_lon"]

# Step 2: Obtain relevant statistics

## Mean duration

mean_duration = stats.mean(duration)
mean_s_lat = stats.mean(s_lat)
mean_s_lon = stats.mean(s_lon)

# Plot data in scatter plots or similar
plt.scatter(s_lon, s_lat, c="r")
plt.scatter(mean_s_lon, mean_s_lat, c="b")
plt.xlabel(r"$\sigma_{lon}$ (deg)")
plt.ylabel(r"$\sigma_{lat}$ (deg)")
plt.savefig("events_statistics.pdf")
