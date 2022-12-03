# In this program we will plot some solar wind properties
# in netCDF format. 

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import netCDF4 as nc
import argparse
import glob

parser = argparse.ArgumentParser(description = "Choose a file to work")
parser.add_argument("--edate", type=str, default="2019-06-22", help="Choose event date")
parser.add_argument("--date", type=str, default="2019-06-22", help="Choose data date")

cmd_args = parser.parse_args()
edate = cmd_args.edate
date = cmd_args.date
year, month, day = date.split("-")

#Load netCDF files

directory = "./solar_data/{}/".format(edate)
format_dir = {"fn1":"f1m", "fn2":"m1m"}
fn1 = directory+"oe_{}_dscovr_s{}{}{}*.nc".format(format_dir["fn1"], year, month, day)
fn2 = directory+"oe_{}_dscovr_s{}{}{}*.nc".format(format_dir["fn2"], year, month, day)
f1 = glob.glob(fn1)[0]
f2 = glob.glob(fn2)[0]
ds1 = nc.Dataset(f1)
ds2 = nc.Dataset(f2)

#Extract relevant data

  # Time column is measured since 1970-01-01 and is measured in miliseconds, so we need to 
  # convert to hours and make sure we first element in our arrays is 0:00:00 
t_raw1 = ds1["time"][:]
t_raw2 = ds2["time"][:]
time1 = 1e-3*(t_raw1 - t_raw1[0])/3600.
time2 = 1e-3*(t_raw2 - t_raw2[0])/3600.
p_speed = ds1["proton_speed"][:]
p_density = ds1["proton_density"][:]
p_temperature = ds1["proton_temperature"][:]
proton_mass = 1.67
units_conversion = 1e-6 # includes ten power of proton mass and conversion factors so the pressure is measured 
                        # in nano pascals
p_dyn_pressure = 0.5*proton_mass*p_density*p_speed**2*units_conversion # dynamic pressure is not included with
                                                                       # date: must be estimated

  # Extract magnetic field components in Geocentric Solar Magnetospheric coordinates
bx = ds2["bx_gsm"][:]
by = ds2["by_gsm"][:]
bz = ds2["bz_gsm"][:]
bt = ds2["bt"][:]

# Make four graphs: Magnetic field versus time, proton density and dynamic pressure versus time,
# proton speed versus time and proton temperature versus time. We will make these graphs separately
# since we may not need all and for having more flexibility when included in paper.

# General settings

sns.set_style("whitegrid")
xticks = [3, 6, 9, 12, 15, 18, 21]
xtick_labels = ["3:00:00", "6:00:00", "9:00:00", "12:00:00", "15:00:00", "18:00:00", "21:00:00"]
xlabel = "UT (hours)"
outdir = "./solar_wind_figures/{}/".format(edate)
  # Magnetic field versus time
print("Plot Magnetic field vs time")
plt.plot(time2, bx, label="bx_gsm")
plt.plot(time2, by, label="by_gsm")
plt.plot(time2, bz, label="bz_gsm")
plt.plot(time2, bt, label="bt")
plt.legend()
plt.xlim(0, 24.)
plt.ylim(-13, 13)
plt.xlabel(xlabel)
plt.ylabel("Mag (nT)")
plt.title("Interplanetary Field Strength Bt and components Bx, By, Bz in Geocentric Solar Magnetospheric coordinates")
plt.xticks(xticks, xtick_labels)
plt.gcf().set_size_inches(10, 4)
plt.tight_layout()
plt.savefig(outdir+"Mfield-{}.pdf".format(date))
plt.clf()

  # Proton speed
print("Plot proton speed ve time")
plt.plot(time1, p_speed, label="proton speed")
plt.legend()
plt.xlim(0, 24.)
plt.ylim(200, 1200)
plt.xticks(xticks, xtick_labels)
plt.xlabel(xlabel)
plt.ylabel(r"Speed $\mathrm{(km~s^{-1})}$")
plt.title("Solar wind proton speed")
plt.tight_layout()
plt.gcf().set_size_inches(10, 3)
plt.savefig(outdir+"Pspeed-{}.pdf".format(date))

  # Proton density and dynamic pressure
print("Plot proton density vs time")
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax1.plot(time1, p_density, label="proton density")
ax2.plot(time1, p_dyn_pressure, "r", label="Approx. dynamic pressure")
ax1.set_xticks(xticks)
ax1.set_xticklabels(xtick_labels)
ax1.set_yticks([10, 20, 30, 40, 50])
ax2.set_yticks([0.0, 1.2, 2.4, 3.6, 4.8, 6.0])
ax1.set_xlim(0., 24.)
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r"Density $\mathrm{(cm^{-3})}$")
ax2.set_ylabel("Dynamic pressure (nPa)")
ax1.set_title("Solar wind proton density and dynamic pressure")
fig.set_size_inches(10, 3)
fig.savefig(outdir+"Pdensity-{}.pdf".format(date))
plt.clf()
  # Proton temperature
print("Print proton temperature vs time")
plt.plot(time1, p_temperature, label="proton temperature")
plt.legend()
plt.xticks(xticks, xtick_labels)
plt.ticklabel_format(axis="y", style="sci")
plt.yscale("log")
plt.xlim(0., 24.)
plt.ylim(1e3, 1e7)
plt.xlabel(xlabel)
plt.ylabel("Temperature (K)")
plt.title("Solar wind proton temperature")
plt.tight_layout()
plt.gcf().set_size_inches(10, 3)
plt.savefig(outdir+"Ptemperature-{}.pdf".format(date))
print("Finished")
