# Mexico map plotter
# The main idea of this program was taken from 
# https://towardsdatascience.com/mapping-with-matplotlib-pandas-geopandas-and-basemap-in-python-d11b57ab5dac
# By Ashwani Dhankhar 
# And the shape file for Mexico from CONABIO
# http://www.conabio.gob.mx/informacion/metadata/gis/destdv250k_2gw.xml?_xsl=/db/meadata/xsl/fgdc_html.xsl&_indent=no

import seaborn as sns
import numpy as np
import pandas as pd
import shapefile as shp
import matplotlib.pyplot as plt
from plotfullmap import plot_map
import argparse
from astropy.table import Table
import glob
import matplotlib.cm as cm

# set figure style
sns.set_style("whitegrid") 
sns.mpl.rc("figure", figsize=(10,6))

# Read shape file of Mexico map
sf = shp.Reader("map.shp")
plot_map(sf)

# Read stations positions table

stations_pos = Table.read("station_data.tab", format="ascii")

# Plot stations positions

plt.plot(stations_pos["Longitude"], stations_pos["Latitude"], "ro")
for i in range(len(stations_pos["Site"])):
    plt.annotate(stations_pos["Site"][i], (stations_pos["Longitude"][i], stations_pos["Latitude"][i]),
		 textcoords="offset points", color="w", xytext=(5, 5), ha="center",
		 bbox=dict(boxstyle="round", pad=0.5, fc="b", alpha=0.7))

ax = plt.gca()
ax.set_aspect("equal", adjustable="box")
plt.savefig("stations_map.pdf")
