import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.table import Table
import pandas as pd

# This program pretends to plot the median Kp index of each event date and the two
# previous days versus date (and use errorbars to estimate variations along these days)
# Uses swarmplot (if i manage to do this kind of plot i will change the energy boxplot for this)

sns.set_style("whitegrid")
table = pd.read_csv("Kp_table.csv")
max_value = table["Max value"]
#colordict = {0:"green", 1:"green", 2:"green", 3:"green", 4:"yellow", 5:"red", 6:"red", 7:"red", 8:"red" , 9:"red"}
color_mask_1 = max_value < 4
color_mask_2 = max_value == 4
color_mask_3 = max_value > 4
#ax sns.swarmplot(x="Date", y="Median", data=table)
x = np.arange(len(table["Date"]))
y = table["Median"]
yerr1 = np.array([table["Variation"][color_mask_1], max_value[color_mask_1]-y[color_mask_1]])
yerr2 = np.array([table["Variation"][color_mask_2], max_value[color_mask_2]-y[color_mask_2]])
yerr3 = np.array([table["Variation"][color_mask_3], max_value[color_mask_3]-y[color_mask_3]])
plt.errorbar(x[color_mask_1], y[color_mask_1], fmt="o", yerr=yerr1, c="green")
plt.errorbar(x[color_mask_2], y[color_mask_2], fmt="o", yerr=yerr2, c="yellow")
plt.errorbar(x[color_mask_3], y[color_mask_3], fmt="o", yerr=yerr3, c="red")
xticks = [table["Date"][0], "", "", "", "", table["Date"][25], "", "", table["Date"][40]]
#xticks[0] = table["Date"][0]
#xticks[-1] = table["Date"][-1]
plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40], labels=xticks)
plt.xlabel("Event date")
plt.ylabel("Kp index median")
plt.title("")
plt.savefig("Kp_index_dist.pdf")
