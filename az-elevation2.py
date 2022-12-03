# This program is a copy of az-elevation.py with the difference that in this case we highlight the satellite positions where TIDs
# were detected in hope of tracing the TIDs origin and propagation

import argparse
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import glob
import pandas as pd
import seaborn as sns

def Change_coords(slon, slat, lon, lat):
    """
    Change coordinate from Latitude-Longitude reference system to local Azimuth-Elevation system
    Inputs:
      slon, slat: Coordinates of station (in degrees)
      lon, lat: Coordinates of object (can be an scalar or array, in degrees)
    Output: AZ, ELE, Coordinates of object in local reference system (in degrees)
    """
    HA = slon - lon
    sin_ELE = np.sin(np.radians(lat))*np.sin(np.radians(slat))+np.cos(np.radians(lat))*np.cos(np.radians(slat))*np.cos(np.radians(HA))
    cos_ELE = np.sqrt(1.-sin_ELE**2)
    cosA = (np.sin(np.radians(lat))-sin_ELE*np.sin(np.radians(slat)))/(cos_ELE*np.cos(np.radians(slat)))
    ELE = np.degrees(np.arcsin(sin_ELE))
    A = np.degrees(np.arccos(cosA))
    AZ = []
    for a, ha in zip(A, HA):
        pi_mask = np.sin(np.radians(ha)) <0
        if pi_mask:
            az = a
        else:
            az = 360. - a
        AZ.append(az)
    return np.array(AZ), ELE

def change_coords2(slon, slat, lon, lat, h):
    """
    Change coordinate from Latitude-Longitude reference system to local Azimuth-Elevation system
    Inputs:
      slon, slat: Coordinates of station (in degrees)
      lon, lat: Coordinates of object (can be an scalar or array, in degrees)
      h: bolide height (in kilometers)
    Output: AZ, ELE, Coordinates of object in local reference system (in degrees)
    In this case we are doing an earth flat aproximation because the equatorial to horizontal transformation works under a wrong assumption
    """
    delta_lat = np.radians(lat - slat)
    delta_lon = np.radians(lon - slon) # delta_lat and delta_lon are defined in a such way that Azimuth is defined 0 degrees at north and growing counterclockwise
    
    Az2 = np.degrees(np.arctan2(delta_lon, delta_lat))
    Re = 6371
    r = Re*np.hypot(delta_lon, delta_lat) # r is the projected distance to fireball measured from station position into the flat-earth plane  
    sin_ele = h /np.hypot(r, h)
    Ele2 = np.degrees(np.arcsin(sin_ele))
    return Az2, Ele2

#  Read parameters from command line and read *.Cmn file
parser = argparse.ArgumentParser(
	  description=""" Choose a file to work""")
parser.add_argument("--date", type=str, default="2000-01-01", help="choose date")
parser.add_argument("--sets", type=str, default="1", help="Choose set of stations")
parser.add_argument("--station", type=str, default="bara",help="Choose station")
parser.add_argument("--Next", action="store_true", help="Use next day data")


cmd_args = parser.parse_args()
date = cmd_args.date
set_folder = cmd_args.sets
Next = cmd_args.Next
station = cmd_args.station.lower()
directory = "./data/{}/set{}/".format(date, set_folder)
#h_satellites_file = "./set{}_TID_detection.csv".format(set_folder)

if Next:
    directory = directory+"next/"
File = glob.glob(directory+"{}*.Cmn".format(station))[0]

f = open(File, "r")
for i in range(4):
    f.readline()
data = f.readlines()
tab = Table.read(data, format="ascii")
sns.set_style("whitegrid")

# Station parameters
s_data = pd.read_csv("./station_data.csv")
s_mask = s_data["Site"] == station.upper()
s_lat = np.array(s_data["Latitude"][s_mask])[0]
s_lon = np.array(s_data["Longitude"][s_mask])[0]

# Meteor parameters
Lat0 = {"2019-06-22":14.9}
Lon0 = {"2019-06-22":-66.2}
vx = {"2019-06-22":-13.4}
vy = {"2019-06-22":6.0}
vz = {"2019-06-22":2.5}
duration = {"2019-06-22":4.5}
h = {"2019-06-22":25}
v_eq = vx[date]
v_lat = vz[date]*np.cos(np.radians(Lat0[date])) - vy[date]*np.sin(np.radians(Lat0[date]))
t = np.linspace(-60, duration[date])
Re = 6371
x,y = Lon0[date] + np.degrees(v_eq/(Re+h[date])*t), Lat0[date] + np.degrees(v_lat/(Re+h[date])*t)


mh = 25 # Meteor height
# Load GLM data
GLM_data = pd.read_csv("./data/{}/GLM16-data.csv".format(date), header=9)
GLM_lat, GLM_lon = GLM_data["latitude"], GLM_data["longitude"]

# Change reference system from latitude-longitude to Azimuth-elevation
Az, Ele = Change_coords(s_lon, s_lat, x, y)
GLM_Az, GLM_ele = Change_coords(s_lon, s_lat, GLM_lon, GLM_lat)

# Now using Earth flat aproximation
az2, ele2 = change_coords2(s_lon, s_lat, x, y, mh)
GLM_az2, GLM_ele2 = change_coords2(s_lon, s_lat, GLM_lon, GLM_lat, h[date])

# Get satellites coordinates and start plotting
AZ, ELE = tab["Az"], tab["Ele"]
# Is useless to plot the whole trajectory of satellites, is only neccesary to plot the trajectory of satellites just a little before fragmentation and after
# Save start time of events
starttime = {"2019-06-22":21 + 25/60.+48/3600.}


time = tab["Time"]

# Load highlighted satellites info

#h_satellites = pd.read_csv(h_satellites_file)
#hstation_mask = h_satellites["Station"] == station.upper()
#h_low = h_satellites["Lower time"][hstation_mask]
#h_up = h_satellites["Upper time"][hstation_mask]


time_window = (time > max(starttime[date]-2.5, 0.0)) & (time < (min(starttime[date]+3.5, 24.0)))
suptitle = "Satellites and meteor trajectory recorded by station {}".format(station.upper())
fig = plt.figure()
prn_list = [2, 5, 6, 9, 12, 13, 17, 19]#np.unique(tab["PRN"])
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]
ax = fig.add_subplot(1,1,1, projection="polar")
outfile = "./az-ele/{}/azimuth-elevation-map-{}-polar-hlght.png".format(date, station)
ax.text(np.radians(0.5*45), 2.1, "Elevation (deg)", fontsize=14)
ax.set_yticks([0, 2./3, 4./3, 2.])
ax.set_yticklabels([r"$90^\circ$", r"$60^\circ$", r"$30^\circ$", r"$0^\circ$"])
tGLM = np.radians(GLM_Az)
rGLM = -2./(0.5*np.pi)*(np.radians(GLM_ele)-0.5*np.pi)
tGLM2 = np.radians(GLM_az2)
rGLM2 = -2./(0.5*np.pi)*(np.radians(GLM_ele2)-0.5*np.pi)
t3 = np.radians(Az)
t2 = np.radians(az2)
r3 = -2./(0.5*np.pi)*(np.radians(Ele)-0.5*np.pi)
R2 = -2./(0.5*np.pi)*(np.radians(ele2)-0.5*np.pi)
if Next: 
    time_window = time < 8.
    suptitle = "Satellites and meteor trajectory recorded by station {} (next day)".format(station.upper())
    outfile = "./az-ele/{}/azimuth-elevation-map-{}-polar-next-hlght.png".format(date, station)
k=0

#h_prn_list = np.unique(h_satellites["PRN"][hstation_mask])
#h_times = h_satellites["Lower time"][hstation_mask]
#h_colors = sns.color_pallete("viridis", len(h_times))
for j in prn_list:
    prn_mask = tab["PRN"] ==j
    theta = np.radians(AZ[prn_mask & time_window])
    r = -2./(0.5*np.pi)*(np.radians(ELE[prn_mask & time_window])-0.5*np.pi)
    if len(theta) > 0:
        ax.plot(theta, r, "o", ms=0.8)#, label="PRN {}".format(j))
        ax.plot(theta[0], r[0], "o", ms=5, c=colors[k])
        ax.text(theta[-1], r[-1], "PRN {}".format(j), c="w", bbox=dict(boxstyle="round", facecolor=colors[k], alpha=0.5))
    else:
        continue
    k+=1
    if k>9:
        k=0 
#for hp in h_prn_list:
#    hl = np.array(h_low[h_satellites["PRN"] ==hp])[0]
#    hu = np.array(h_up[h_satellites["PRN"] ==hp])[0]
#    h_time_window = (time <hu) & (time >hl)
#    h_prn_mask = tab["PRN"] == hp
#    h_prn_mask2 = h_satellites["PRN"] == hp
#    theta2 = np.radians(AZ[h_prn_mask & h_time_window])
#    r2 = -2./(0.5*np.pi)*(np.radians(ELE[h_prn_mask & h_time_window])-0.5*np.pi)
#    if len(theta2) > 0:
#        ax.scatter(theta2[0], r2[0], c="y")#, ms=10)
#    else:
#        continue
#fig.colorbar()
#ax.plot(tGLM, rGLM, "m.", label="GLM data")
ax.plot(tGLM2, rGLM2, "k.", label="GLM data")
#ax.plot(t3, r3, "k--", label="Meteor trajectory")
#print(len(t2), len(R2))
ax.plot(t2, R2, "m--", label="Meteor trajectory")
legend =ax.legend()
legend.get_frame().set_alpha(0.5)
ax.set_xlabel("Azimuth (deg)", fontsize=14)
fig.suptitle(suptitle)
fig.set_size_inches(8,6)
plt.savefig(outfile)
