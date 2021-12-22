import numpy as np
from scipy.interpolate import griddata
import pandas as pd
import matplotlib.pyplot as plt
import shapefile as shp
from plotfullmap import plot_map
import matplotlib.cm as cm
import argparse

parser = argparse.ArgumentParser(
    description=""" Choose a file to work""")


parser.add_argument('--date', type=str, default='2000-01-01',
			       help='Choose date. Format: yyyy-mm-dd')

parser.add_argument("--previous", action="store_true", help="Choose previous days data")

cmd_args = parser.parse_args()
date = cmd_args.date
previous = cmd_args.previous
directory = "./data/"+date
if previous:
    directory = "./data/"+date+"/previous"
df = pd.read_csv(directory+"/detrended_TEC.csv")
latitude = df["Latitude"]
longitude = df["Longitude"]
det_TEC = df["vTEC"]

xmin, ymin = min(longitude), min(latitude)
xmax, ymax = max(longitude), max(latitude)
res = 0.2
Nx, Ny = (xmax-xmin)/res, (ymax-ymin)/res
grid_x, grid_y = np.mgrid[xmin:xmax:Nx*1j, ymin:ymax:Ny*1j]

fig =plt.figure()

# Read and plot shape file of Mexico map
sf = shp.Reader("map.shp")
plot_map(sf)

#Plot GLM trajectories

data_16 = pd.read_csv(directory+"/GLM/GLM-16-data.csv", header=9)
data_17 = pd.read_csv(directory+"/GLM/GLM-17-data.csv", header=9)

plt.plot(data_16["longitude"], data_16["latitude"], "r.")
plt.plot(data_17["longitude"], data_17["latitude"], "r.")


# Plot event position

event = pd.read_csv("GLM-database.csv")
mask = event["Fecha"] == date
event_lat = event["Latitud"][mask]
event_lon = event["Longitud"][mask]
event_ID = str(event["ID"][mask])
#event_lat, event_lon =  22.5, -83.8
plt.plot(event_lon, event_lat, "m*")
plt.annotate(event_ID, (event_lon, event_lat), textcoords="offset points", color="w", fontsize="small",
    xytext=(10, 10), ha="center", bbox=dict(boxstyle="round", pad=0.5, fc="m", alpha=0.7))

#vx, vy = -2.4, 13.6
#height = 23.7 


points = []
for x, y in zip(longitude, latitude):
    points.append([x, y])

#grid_z0 = griddata(points, det_TEC, (grid_x, grid_y), method="nearest")
#grid_z1 = griddata(points, det_TEC, (grid_x, grid_y), method="linear")
grid_z2 = griddata(points, det_TEC, (grid_x, grid_y), method="cubic")

#plt.subplot(221)
#plt.imshow(grid_z0.T, extent=(xmin, xmax, ymin, ymax), origin="lower")

#plt.subplot(222)
#plt.imshow(grid_z1.T, extent=(xmin, xmax, ymin, ymax), origin="lower")

#plt.subplot(223)
im = plt.imshow(grid_z2.T, extent=(xmin, xmax, ymin, ymax), origin="lower")

cbar = plt.colorbar(im)
cbar.set_label("Delta vTEC (TECU)")
plt.title("Detrended TEC for 1 station for 2019-02-01")
plt.gcf().set_size_inches(6, 6)
outdir = "./TEC_grid_maps/"
plt.savefig(outdir+date+"_TEC_map.pdf")
