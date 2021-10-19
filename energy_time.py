import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import seaborn as sns

# This program pretends to plot The kinetic energy of the meteors of the sample
# versus the respecting date.

# Read meteors database table

t = Table.read("meteors_database.tab", format="ascii")
t2 = Table.read("USG_meteors_database.tab", format="ascii")
# extract relevant data

ID , date, energy = t["ID"], t["Fecha"], t["Total energy recallibrated (kT)"]
ID_USG, date_USG, energy_USG = t2["ID"], t2["Fecha"], t2["Energia total (kT)"]
energy_arr = []
de = []

for e in energy:
    a, b = e.split("+/-")
    energy_arr.append(float(a))
    de.append(float(b))


# Make Boxplot (with whiskers)
sns.set_style("whitegrid")
GLM_energy = np.array(energy_arr)
USG_energy = t2["Energia total (kT)"]

plt.boxplot([GLM_energy, USG_energy], showfliers=False, patch_artist=True, 
             boxprops=dict(facecolor="red", color="red", alpha=0.5),         
             whiskerprops=dict(color="red"),
             capprops=dict(color="red"), 
             medianprops=dict(color="blue"))
plt.xticks([1, 2], ["GLM", "USG"])
plt.xlabel("Meteors sample")
plt.ylabel("Total energy (kT)")
plt.title("Meteors energy distribution (outliers removed)")
plt.savefig("energies_boxplot.pdf")
