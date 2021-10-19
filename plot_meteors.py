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

# plot positions plus uncertainties in the map

plt.errorbar(f["Longitud"], f["Latitud"], xerr=f["sig_lon"], yerr=f["sig_latitude"], fmt="bo", capsize=3)
# Offset of labels
x_off = [10, 10, 10, -10, -10, -10, 10, 10, 10, 10, 10, -10, -10, 10, 10, 10, 10, 10, 10, 10, -10, -10, 10, 10, 10, -10, 10, -10, -10, 10, -10, 10, -10, 10, 10]
y_off = [10, 10, 10, -10, 10, -10, 10, 10, 10, 10, -10, -10, -10, 10, 10, 10, -10, -10, 10, -10, 10, 10, 10, 10, 10, -10, -10, 10, -10, 10, 10, 10, 10, 10, 10]
for i in range(len(f["ID"])):
    plt.annotate(f["ID"][i], (f["Longitud"][i], f["Latitud"][i]), textcoords="offset points", color="w", fontsize="small",
    xytext=(x_off[i], y_off[i]), ha="center", bbox=dict(boxstyle="round", pad=0.5, fc="b", alpha=0.7))
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.savefig("meteors_map.pdf")
