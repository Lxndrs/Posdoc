import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# The goal of this program is to estime the discrepancy factor between GLM energy estimations
# and the energy estimations by USG instruments in order to calibrate uncertainties to other
# events only registered by the GLM. For this we will make a linear fit between GLM and USG
# energies for events where data is available in both the GLM website and CNEOS fireballs website
# (Only three events meet this criteria: GLM-00, GLM-23 and GLM-Ven).
# The slope will be the correction factor and the residuals the uncertainties. 

sns.set_style("whitegrid")

# insert meteor energies

GLM_energies = np.array([1.0763569, 0.065099193, 3.1305469])
USG_energies = np.array([1.4, 0.076, 6])

GLM_uncertainties = np.array([0.18966796, 0.026813257, 0.28006880])
weights = 1./(1+GLM_uncertainties)

# Do the fit: p is an array with the fit coefficients, res is the residuals, rank, single values 
# and rcond are other outputs not used
p, res, rank, s_values, rcond = np.polyfit(GLM_energies, USG_energies, 1, w=weights, full=True)

# make a linar polynomial with the coefficients
poly = np.poly1d(p)
#poly_plus = np.poly1d([p[0]+res[0], p[1]])
#poly_minus = np.poly1d([p[0]-res[0], p[1]])

# make the plot
f = plt.figure()
ax1 = f.add_subplot(1, 2, 1)
ax2 = f.add_subplot(1, 2, 2)
ax1.errorbar(GLM_energies, USG_energies, xerr=GLM_uncertainties, fmt="g.", label="data points")
ax1.plot(GLM_energies, poly(GLM_energies), "k-", label="linear fit")
#ax1.set_xlabel("Meteor estimated kinetic energy by GLM (kilotons)", fontsize="small")
ax1.set_ylabel("Meteor estimated kinetic energy by USG sensors (kilotons)", fontsize="small")
ax1.text(0, 6.5, r"Linear fit: $y={:.3f}x - {:.3f}$".format(p[0], np.abs(p[1])), fontsize="small")
ax2.plot(GLM_energies, USG_energies-poly(GLM_energies), "go")
ax2.plot([GLM_energies[0], GLM_energies[0]], [0, USG_energies[0]-poly(GLM_energies[0])], "b--")
ax2.plot([GLM_energies[1], GLM_energies[1]], [0, USG_energies[1]-poly(GLM_energies[1])], "b--")
ax2.plot([GLM_energies[2], GLM_energies[2]], [0, USG_energies[2]-poly(GLM_energies[2])], "b--")
ax2.text(0, 1.05, "residuals (squares sum) = {:.4f}".format(res[0]), fontsize="small")
ax2.axhline(y=0, color="k")
ax2.set_ylim(-1, 1)
#ax2.set_xlabel("Meteor estimated kinetic energy by GLM (kilotons)", fontsize="small")
ax2.set_ylabel("Residuals")
f.text(0.5, 0.01, "Meteor estimated kinetic energy by GLM (kilotons)", ha='center', fontsize="small")
#plt.fill_between(GLM_energies, poly_minus(GLM_energies), poly_plus(GLM_energies), color="r", alpha = 0.2)
#plt.plot(GLM_energies, poly_minus(GLM_energies), "r--", alpha=0.2)
#plt.plot(GLM_energies, poly_plus(GLM_energies), "r--", alpha=0.2)
#plt.legend()
f.tight_layout()
plt.savefig("energy_fit.pdf")
