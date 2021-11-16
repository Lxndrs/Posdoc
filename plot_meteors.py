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


# Read meteors_database

f = Table.read("meteors_database.tab", format="ascii")
f2 = Table.read("USG_meteors_database.tab", format="ascii")
# plot positions plus uncertainties in the map
GLM_mask = f["ID"] == "GLM-Ven" # We will exclude the Venezolan meteor since it will distort the map
USG_mask = f2["ID"] == "USG-09" # We will exclude the Venezolan meteor since it will distort the map

plt.errorbar(f["Longitud"][~GLM_mask], f["Latitud"][~GLM_mask], xerr=f["sig_lon"][~GLM_mask], yerr=f["sig_latitude"][~GLM_mask], fmt="bo", capsize=3)
plt.plot(f2["Longitud"][~USG_mask], f2["Latitud"][~USG_mask], "go")
# Offset of labels
x_off = [10, 10, 10, -10, -10, -10, 10, 10, -10, 10, 10, -10, -10, 10, 10, 10, 10, 10, 10, 10, -10, -10, 10, 10, 10, -10, 10, -10, -10, 10, -10, 10, -10, 10, 10]
y_off = [10, 10, 10, -10, 10, -10, 10, 10, -10, 10, -10, -10, -10, 10, 10, 10, -10, -10, 10, -10, 10, 10, 10, 10, 10, -10, -10, 10, -10, 10, 10, 10, 10, 10, 10]

x2_off = [10, -10, 10, 10, 10, -10, 10, -10, 10]
y2_off = [10, 10, -10, -10, 10, -10, 10, -10, -10]

for i in range(len(f["ID"][~GLM_mask])):
    plt.annotate(f["ID"][~GLM_mask][i], (f["Longitud"][~GLM_mask][i], f["Latitud"][~GLM_mask][i]), textcoords="offset points", color="w", fontsize="small",
    xytext=(x_off[i], y_off[i]), ha="center", bbox=dict(boxstyle="round", pad=0.5, fc="b", alpha=0.7))
for i in range(len(f2["ID"][~USG_mask])):
    plt.annotate(f2["ID"][~USG_mask][i], (f2["Longitud"][~USG_mask][i], f2["Latitud"][~USG_mask][i]), textcoords="offset points", color="w", fontsize="small",
    xytext=(x2_off[i], y2_off[i]), ha="center", bbox=dict(boxstyle="round", pad=0.5, fc="g", alpha=0.7))

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.savefig("meteors_map.pdf")
