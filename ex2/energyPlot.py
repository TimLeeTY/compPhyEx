"""
=============================
Plotting the change in Energy
=============================
The θ and ⍵ values are loaded from the file 'energyDiffOut.csv' which lets us
find the corresponding energies
We then find the fractional change in energy which we plot on separate axis of
the same figure due to the vast difference in scales
"""
import matplotlib.pyplot as plt
import numpy as np


def E(y): return(-1*np.cos(y[0])+y[1]**2/2)  # finds energy based on θ and ⍵


# importing variables
tarr, thetaRK, omegaRK, thetaLS, omegaLS = np.loadtxt(
    "energyDiffOut.csv", unpack=True, delimiter=',')
E_RK, E_LS = E(np.array([thetaRK, omegaRK])), E(np.array([thetaLS, omegaLS]))

fig, ax1 = plt.subplots()
plt.xlim(0, tarr[-1]*10**(-3)/(2*np.pi))  # scale x axis by 10^4
plt.xlabel(r"$t/T_\mathrm{0}\times 10^{3}$")
ax1.set_ylabel(
    r'$\Delta E/E_0\times 10^{-9}$ (4th order R-K)', color='#1f77b4')
ax1.set_ylim(np.array([0, 14]))
ax1.tick_params('y', colors='#1f77b4')
# plot of total energy against time (relative to starting energy)
ax1.plot(tarr*10**(-3)/(2*np.pi),
         (E_RK-E_RK[0])/E_RK[0]*10**(9), label="4th order R-K")
ax2 = ax1.twinx()  # plot either result on separate axes of same graph (they differ by ~10^3)
ax2.plot(tarr*10**(-3)/(2*np.pi),
         (E_LS-E_LS[0])/E_LS[0]*10**(6), label='LSODA', color='#ff7f0e')
ax2.set_ylabel(r'$\Delta E/E_0\times 10^{-6}$ (LSODA)', color='#ff7f0e')
ax2.set_ylim(np.array([0, 10]))
ax2.tick_params('y', colors='#ff7f0e')
fig.savefig("energyDiff.pdf", format="pdf")
plt.show()
