# The purpose of this program is to evaluate the azimuth-elevation plots for all satellites which could make contact with
#our GPS stations and took data of TEC. Also we now seek to constrain the TIDs velocities in order to find a a possible range

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
#parser.add_argument("--polar", action="store_true", help="Plot maps in polar coordinates")

cmd_args = parser.parse_args()
date = cmd_args.date
set_folder = cmd_args.sets
Next = cmd_args.Next
station = cmd_args.station.lower()
directory = "./data/{}/set{}/".format(date, set_folder)
#polar = cmd_args.polar
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
t = np.linspace(0, duration[date])
Re = 6371
x,y = Lon0[date] + np.degrees(v_eq/(Re+h[date])*t), Lat0[date] + np.degrees(v_lat/(Re+h[date])*t)

# Plot propagation circle
v_TID_high = 0.362
delta_t = 3600
r_circ = v_TID_high*delta_t
r_circ2 = v_TID_high*2*delta_t
if Next:
    r_circ = r_circ*3
    r_circ2 = r_circ2*2
Re = 6371

R = np.degrees(r_circ/(Re+h[date]))
R2 = np.degrees(r_circ2/(Re+h[date]))
x_circ, y_circ = R*np.cos(np.linspace(-np.pi, np.pi))+Lon0[date], R*np.sin(np.linspace(-np.pi, np.pi))+Lat0[date]
x_circ2, y_circ2 = R2*np.cos(np.linspace(-np.pi, np.pi))+Lon0[date], R2*np.sin(np.linspace(-np.pi, np.pi))+Lat0[date]
# Load GLM data
GLM_data = pd.read_csv("./data/{}/GLM16-data.csv".format(date), header=9)
GLM_lat, GLM_lon = GLM_data["latitude"], GLM_data["longitude"]

# Change reference system from latitude-longitude to Azimuth-elevation
Az, Ele = Change_coords(s_lon, s_lat, x, y)
Az_flat, Ele_flat = change_coords2(s_lon, s_lat, x, y, h[date])
GLM_Az, GLM_ele = Change_coords(s_lon, s_lat, GLM_lon, GLM_lat)
GLM_Az2, GLM_ele2 =change_coords2(s_lon, s_lat, GLM_lon, GLM_lat, h[date])
az_circ, ele_circ = change_coords2(s_lon, s_lat, x_circ, y_circ, h[date])
az_circ2, ele_circ2 = change_coords2(s_lon, s_lat, x_circ2, y_circ2, h[date])

# Get satellites coordinates and start plotting
AZ, ELE = tab["Az"], tab["Ele"]
# Is useless to plot the whole trajectory of satellites, is only neccesary to plot the trajectory of satellites just a little before fragmentation and after
# Save start time of events
starttime = {"2019-06-22":21 + 25/60.+48/3600.}

# Estimate the time that satellites could have recorded data (before the 24.0 hrs and thus available in the 
# original data set)
remaining_time = 24.0 - starttime[date]
# End day circle
end_day_circ = v_TID_high*remaining_time*delta_t
end_day_radius = np.degrees(end_day_circ/(Re+h[date]))
x_end, y_end = end_day_radius*np.cos(np.linspace(-np.pi, np.pi))+Lon0[date], end_day_radius*np.sin(np.linspace(-np.pi, np.pi))+Lat0[date]
az_circ_end, ele_circ_end = change_coords2(s_lon, s_lat, x_end, y_end, h[date])

time = tab["Time"]
time_window = (time > max(starttime[date]-1.1, 0.0)) & (time < (min(starttime[date]+3.5, 24.0)))
label1 = "1 hour TID propagation"
label2= "2 hours TID propagation"
suptitle = "Satellites and meteor trajectory recorded by station {}".format(station.upper())
outfile = "./az-ele/{}/azimuth-elevation-map-{}.pdf".format(date, station)
fig = plt.figure()
prn_list = np.unique(tab["PRN"])
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]

ax = fig.add_subplot(1,1,1, projection="polar")
outfile = "./az-ele/{}/azimuth-elevation-map-{}-polar.pdf".format(date, station)
#    ylabel = ""
ax.text(np.radians(0.5*45), 2.1, "Elevation (deg)", fontsize=14)
ax.set_yticks([0, 2./3, 4./3, 2.])
ax.set_yticklabels([r"$90^\circ$", r"$60^\circ$", r"$30^\circ$", r"$0^\circ$"])
t1 = np.radians(az_circ)
r1 = -2./(0.5*np.pi)*(np.radians(ele_circ)-0.5*np.pi)
t2 = np.radians(az_circ2)
r2 = -2./(0.5*np.pi)*(np.radians(ele_circ2)-0.5*np.pi)
t_end = np.radians(az_circ_end)
r_end = -2./(0.5*np.pi)*(np.radians(ele_circ_end)-0.5*np.pi)
tGLM = np.radians(GLM_Az)
rGLM = -2./(0.5*np.pi)*(np.radians(GLM_ele)-0.5*np.pi)
tGLM2 = np.radians(GLM_Az2)
rGLM2 = -2./(0.5*np.pi)*(np.radians(GLM_ele2)-0.5*np.pi)
t3 = np.radians(Az)
r3 = -2./(0.5*np.pi)*(np.radians(Ele)-0.5*np.pi)
t_flat, r_flat = np.radians(Az_flat), -2./(0.5*np.pi)*(np.radians(Ele_flat)-0.5*np.pi)

if Next: 
    label1 = "3 hours TID propagation"
    label2= "4 hours TID propagation"
    time_window = time < 8.
    suptitle = "Satellites and meteor trajectory recorded by station {} (next day)".format(station.upper())
    outfile = "./az-ele/{}/azimuth-elevation-map-{}-polar-next.pdf".format(date, station)
#ax.text(t1[-1], r1[-1], label1, c="r")
#ax.text(t2[-1], r2[-1], label2, c="b")
#ax.text(t_end[-1], r_end[-1], "TID propagation at 24 hrs UT", c="g")

k=0
for j in prn_list:
    prn_mask = tab["PRN"] ==j
    theta = np.radians(AZ[prn_mask & time_window])
    r = -2./(0.5*np.pi)*(np.radians(ELE[prn_mask & time_window])-0.5*np.pi)
    if len(theta) > 0:
        ax.plot(theta, r, "o", ms=0.8, label="PRN {}".format(j))
        ax.text(theta[-1], r[-1], "PRN {}".format(j), c=colors[k])
    else:
        continue
    k+=1
    if k>9:
        k=0 
#ax.plot(t1, r1, "r--", label=label1)
#ax.plot(t2, r2, "b--", label=label2)
ax.plot(t_flat, r_flat, "m--")
ax.plot(tGLM, rGLM, "m.", label="GLM data")
ax.plot(tGLM2, rGLM2, "k.", label="GLM data")
ax.plot(Az_flat, Ele_flat, "m--")
ax.plot(t3, r3, "k--", label="Meteor trajectory")
#ax.plot(t_end, r_end, "g--")
ax.set_xlabel("Azimuth (deg)", fontsize=14)
fig.suptitle(suptitle)
fig.set_size_inches(14,12)
plt.savefig(outfile)
