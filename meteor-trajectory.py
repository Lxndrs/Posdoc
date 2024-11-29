# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# By Alejandro Tarango Yong and Raul Gutierrez-Zalapa
#
# The second part in the creation of the caribbean sTEC maps is the estimation of the meteor
# trajectory. We will work with the meteor data available at https://cneos.jpl.nasa.gov/fireballs/
# for caribbean meteor. At thos database, we found the following parameters:
# 
# - Date and time of fragmentation: 2019-06-22, 21:25:48 UT
# - Latitude and longitude: 14.9, -66.2
# - Altitude: 25 km
# - Total velocity: 14.9 km/s
# - Velocity components: v_x = -13.4, v_y = 6.0, v_z = 2.5
#   (Where v_x is parallel to the equatorial plane and towards the prime meridian, v_z is parallel 
#    to the Earth's rotation axis and v_y completes the orthonormal set.)
# - Total radiated energy: 2.8e10 J
# - Calculated total impact energy: 6 kT
#
# The meteor is descending towards Earth's surface with a velocity v(t(h)). The angle of the 
# trajectory respect to the surface is \theta(t(h)). In order to plot the trajectory we will 
# need the meteor position (latitude and longitude) as a function of time.
# First of all, we will need to split the velocity into components:
#
#         (1)  v_h = v\sin\theta       velocity normal to the Earth's surface
#         (2)  v_{tan} = v\cos\theta   velocity parallel to Earth's surface
#
# Only v_{tan} is useful for the map creation. Also, we can split v_{tan} into components: the
# first one in the longitudinal direction and the second one in the latitudinal direction.
# v_{tan} should not change in direction, since the Coriolis force is negligible in this case, then 
# the angle between the trajectory and the local parallel, called \alpha, can be estimated using 
# the velocity components found at any time (including the given by the CNEOS database):
#
#         (3) \tan\alpha = v_{lat0}/v_{lon0}
#         (4) v_{lat0} = v_z\cos(L0) - v_y\sin(L0)  where L0 = 14.9
#         (5) v_lon0 = v_x
#         (6) v_{lat} = v_{tan}\sin\alpha = v\cos\theta\sin\alpha
#         (7) v_{lon} = v_{tan}\cos\alpha = v\cos\theta\cos\alpha
#
# Now, to convert these velocities into position angles, we need to divide by the circumference
# radius (a.k.a the actual height plus the radius of Earth), and since the velocity and height
# are variable, the latitude and position can be estimated by the following **integrals**:

#         (8) dL = (v_{lat}/(R_E + h)) dt
#         (9) d\lambda = (v_{lon}/(R_E + h)) dt
#         
# Now, time and height follow a relation like dt = m dh, where m is the conversion factor, which 
# could be contant or variable, but we found from Raul's work that in fact is a constant. The same
# happens to the relation between h and \theta, so we can do the next:
#
#       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      | (10)              L - L0 = m\sin\alpha\int_{h0}^{h} (v\cos\theta)/(R_E + h) dh |  
#      | (11) \lambda - \lambda_0 = m\cos\alpha\int_{h0}^{h} (v\cos\theta)/(R_E + h) dh | 
#       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# The integral can be solved numerically since we count with data of v and \theta as a function of 
# h for Raul's work, as well as the numerical valus of m. Also, we can extrapolate back in time 
# the meteor trajectory since before the estimations made by Raul the meteor velocity is constant
# and equal to the terminal velocity which can be estimated from Raul's work as well.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Import python libriaries
#

import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table
from scipy.integrate import quad

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Extract data from Raul's work. He only gave me the graphs in form of image files, so I had to
# estimate these vales by visual inspection.
# height units ate kilometers, time units are seconds, velocity units are km/s

h_sample = np.array([70, 60, 50, 40, 30, 25])         # Raul's data started in time when h= 70 km
t_sample = np.array([0, 1.15, 2.3, 3.5, 4.8, 5.5]) 
v_sample = np.array([20.6, 20.6, 20.5, 19.5, 16.5, 15])
theta = np.radians([24.6, 24.82, 25.04, 25.26, 25.5, 25.62]) # Original data for theta is 
                                                             #in degrees

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Incorporate data from CNEOS and obtain new parameters
# 

R_E = 6371 # Radius of Earth (in kilometers)
h0 = 25    # height (in kilometers)
lon0, lat0 = -66.2, 14.9
v_x, v_y, v_z = -13.4, 6.0, 2.5 # velocity components (in km/s), minus in v_x implies to west
v_eq0 = v_x                     # and positive in v_z and v_y implies to north and falling to Earth
v_lat0 = v_z*np.cos(np.radians(lat0)) - v_y*np.sin(np.radians(lat0))
alpha = np.arctan2(v_lat0, v_eq0)
coef1 = np.polyfit(h_sample, t_sample, 1) # Linear fit of time as a function of height
coef2 = np.polyfit(t_sample, h_sample, 1) # The inverse relation
m = coef1[0]                              # We want the slope of the fit

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Interpolation of data. We sould create functions for the integrand in equations (10) and (11) in
# the scenarios when h < 70 km (use Raul's data) and h > 70 km (the mateors moves at constant 
# velocity and theta varies in linear way with height (or time)).

f1 = interp1d(h_sample,v_sample*np.cos(np.radians(theta))/(R_E + h_sample)) # integrand of 
                                                                            # equations (10) & (11)
I1 = np.array([quad(f1, h0, h)[0] for h in h_sample])                       # integral of equations 
                                                                            #(10) and (11)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to compute the actual position of the meteor as a function of time
# The idea is to plot the meteor trajectory until the current time. If the required time is
# before meteor passage, then plot nothing, and if is after fragmentation, thenplot the whole 
# trajectory. Also we must show important events if they occured at the required time (i.e passage
# through ionosphere, fragmetation and detection, at least).

# important events (so far):
# - Passage through the IPP
# - Breakup
# - USG detection

h_IPP = 350    # Ionospheric Piercing Point height
h_breakup = 50 # Breakup height
h_USG = h0
t_vs_h = np.poly1d(coef1)  # Linear realtion between time and height
h_vs_t = np.poly1d(coef2) # The inverse relation (will be very useful)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Set the time intervals --> Convert simulation time into UT. We know that the meteor was deteded 
# at 21:25:48 UT (21.43 UT) and this happened at a height of 25 km. With, this, we can create a 
# time array beween the time when h = 400 km (the ionospheric Piercing point) and h = 25 km, when 
# the USG sensors detected the meteor, get the correspondent heights and calculate the positions 
# with equations (10) and (11). For the terminal-velocty part, we found that the meteor had a 
# terminal velocity of 20.6 km/s
# The conversion between simulation time and UT is
#
#         (12) UT = (t - 5.5)/3600 + 21.43

UT_USG = 21.43
v_ter = 20.6

# Longitude and latitude of meteor near breakup (original arrays), uppercase works for longitude, 
# lowercase for latitude

L1, l1 = lon0 + np.degrees(m*np.cos(alpha)*I1), lat0 + np.degrees(m*np.sin(alpha)*I1) 


# Compute extrapolation for theta(h) and terminal-velocity stage of trajectory
th_coef = np.polyfit(h_sample, theta, 1)
th_vs_h = np.poly1d(th_coef)
h_extra = np.array([400, 350, 300, 250, 200, 150, 100, 80, 70])
t_extra = t_vs_h(h_extra)
h1 = 70
f2 = interp1d(h_extra, np.cos(th_vs_h(h_extra))/(R_E + h_extra))
I2 = np.array([quad(f2, h1, h)[0] for h in h_extra])
L2, l2 = lon0 + np.degrees(m*v_ter*np.cos(alpha)*I2), lat0 + np.degrees(m*v_ter*np.sin(alpha)*I2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# concatenate l1 and l2, L1 and L2
h_out, t_out = np.concatenate((h_extra[:-1], h_sample), axis=None), np.concatenate((t_extra[:-1], 
                               t_sample), axis=None)
L, l = np.concatenate((L2[:-1], L1), axis=None), np.concatenate((l2[:-1], l1), axis=None)

# Compute Universal time array
UT_out = (t_out - 5.5)/3600. + UT_USG

# Compute important event times and positions
t_IPP = (t_vs_h(h_IPP) - 5.5)/3600. + UT_USG
t_breakup = (t_vs_h(h_breakup) - 5.5)/3600. + UT_USG
t_USG = (t_vs_h(h_USG) - 5.5)/3600. + UT_USG
L_int = interp1d(UT_out, L)
l_int = interp1d(UT_out, l)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Write output
#
output_table = Table([UT_out, h_out, L, l], names=("UT", "Height (km)", "Longitude (deg)", 
                                                   "Latitude (deg)"))
outfile = "meteor-trajectory.tab"
output_table.write(outfile, format="ascii", overwrite=True)
print("Meteor trajectory traced succesfully")
print("Time of IPP passage (UT): {}".format(t_IPP))
print("Time of breakup (UT): {}".format(t_breakup))
print("Time of USG detection (UT): {}".format(t_USG))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
