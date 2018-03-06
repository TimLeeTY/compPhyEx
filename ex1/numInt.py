
"""
===================================================
Ex.1 Numerical Integration using Monte-Carlo method
===================================================
"""

import numpy as np
import matplotlib.pyplot as plt


def MCint(N, n):  # N holds number of samples and n is the number of trials per N
    D, s = 8, np.pi/8
    V = s**8
    fmean, fsqmean = np.zeros(n), np.zeros(n)
    for i in range(n):
        # sums 8 uniformly distributed random numbers in [0,s]
        r = (np.random.rand(D, N)).sum(axis=0)*s
        f = np.sin(r)*V*10**6  # evaluates sin(x0+x1+...+x7)
        fmean[i] = f.mean()  # evaluates <f>
        fsqmean[i] = (f**2).mean()  # evaluates <f^2>
    sigma = np.sqrt((fsqmean-fmean**2)/N)  # finds error in estimate
    # return the best estimate of the integral in n trials, the standard deviation
    # of that mean, and the mean of the error estimate
    return(fmean.mean(), fmean.std(), sigma.mean())


#%%
"""
-----------------------
§1. Evaluating integral
-----------------------
Repeats the integral n times for each value of N samples.
Take mean of n trials as the best estimate for the integral.
"""
n = 25
vMC = np.vectorize(MCint)
target = 10**6*(70-16*np.sin(np.pi/8)+56*np.sin(np.pi/4)-112*np.sin(3*np.pi/8))

sample = (2**np.arange(5, 22))
lsamp = np.log(sample)
[y, yerr, zerr] = vMC(sample, n)

plt.figure()
plt.semilogx()
plt.ylim(525, 545)
plt.yticks(np.arange(525, 551, 5))
plt.xlabel(r'$N$')
plt.ylabel(r'Integral value')
plt.errorbar(sample, y, yerr, fmt='+', color='#1f77b4', lw=0.8, capsize=3)
plt.plot(sample, np.ones(len(sample))*target, '-', color='#2ca02c', lw=0.8)

plt.savefig("est.pdf", format="pdf")
#%%
"""
---------------------
§2. Evaluating errors
---------------------
Comparing the error obtained from the standard deviation of n trials to the
average error from equation (3) of the handout
"""
fig, ax = plt.subplots()
ax.loglog(basex=2)
plt.ylim(10**(-5), 10**(-1))
yfit = np.polyfit(lsamp, np.log(yerr/target), 1)
yffit = np.poly1d(yfit)
ax.plot(sample, yerr/target, '+', color='#1f77b4', label=r'S.D. in mean')
ax.plot(sample, np.exp(yffit(lsamp)), '-', color='#1f77b4', linewidth=0.8)

zfit = np.polyfit(lsamp, np.log(zerr/target), 1)
zffit = np.poly1d(zfit)
ax.plot(sample, zerr/target, 'x', color='#ff7f0e', label=r'$\sigma (N)$')
ax.plot(sample, np.exp(zffit(lsamp)), '-', color='#ff7f0e', linewidth=0.8)

plt.xlabel(r'$N$')
plt.ylabel(r'Fractional Error, $\Delta$')
plt.grid(True)
ax.legend()
fig.savefig("err.pdf", format="pdf")
