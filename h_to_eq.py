# This program is aimed to change from horizontal coordinates to equatorial coordinates
# Horizontal coordinates of events are got from AMS data and equatorial coordinates would give us 
# the place where the event could be seen at the zenith. If we get an event with two or more
# witnesses we could find a parallax and thus estimate the altitude

import numpy as np
from scipy.optimize import fsolve
import argparse


def equations(q, Al, Az, Lat):
    """
    Set of equations to solve
    """

    x, y, z = q
    k1 = np.sin(Az)*np.cos(Al)
    k2 = np.cos(Az)*np.cos(Al)*np.sin(Lat) + np.sin(Al)*np.cos(Lat)
    k3 = -np.cos(Az)*np.cos(Al)*np.cos(Lat) + np.sin(Al)*np.sin(Lat)
    return (x*z - k1, y*z - k2, 1 - z**2 - k3**2) 

# Input data

parser = argparse.ArgumentParser(description= "Insert local coordinates of events")
parser.add_argument("--latitude", type=float, help="Local latitude of observer (in degrees)")
parser.add_argument("--longitude", type=float, help="Local longitude of observer (in degrees)")
parser.add_argument("--azimuth", type=float, help="Azimuth of the event (in degrees)")
parser.add_argument("--altitude", type=float, help="Altitude of the event (in degrees)")
parser.add_argument("--h0", type=float, default=0, help="Initial guess of hour angle (in degrees)")
parser.add_argument("-d0", type=float, help="Initial guess of declination (in degrees)")

cmd_args = parser.parser_args()
l_latitude = cmd_args.latitude
azimuth = cmd_args.azimuth
altitude = cmd_args.altitude
h0 = cmd_args.h0
delta0 = cmd_args.d0
l_longitude = cmd_args.longitude

# set initial guess

x0 = np.sin(np.radians(h0))
y0 = np.cos(np.radians(h0))
z0 = np.cos(np.radians(delta0))

# solve the system

sol = fsolve(equations, (x0, y0, z0), args=(np.radians(altitude), np.radians(azimuth), np.radians(l_latitude)))

# translate solution into (latitude, longitude)

latitude = np.degrees(np.arccos(z) -0.5*np.pi) # Arc cosine range is from 0 tp \pi and latitude is shifted 90 degrees

if x>=0: # if sin(h) is positive, then 0 <= h <= 180
    h = np.arccos(y)
else: # if sin(h) is negative, then 180 < h < 360
    h = np.arccos(y) + np.pi


# move to a longitude where h=0 (and thus we may see the event at the zenith)

longitude = l_longitude - np.degrees(h)
