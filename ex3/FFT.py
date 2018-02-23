"""
=====================================================
Ex.3A Apply FFT method to obtain diffraction patterns
=====================================================
"""

import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
# units of microns used consistently

# computes the aperture based on the slit width, number of slits and their separation
# x coordinate scaled to match dimensions of Fourier transform by Delta=L/N


def aper(d, N, L, funcA=lambda x, y, z: 1, slits=1, sep=None):
    if sep == None:
        sep = d
    elif sep < d or slits*sep > L or type(slits) != int:
        print('slit dimensions invalid')
        return()
    Delta = (L/N)
    x = np.arange(N)*Delta
    fA = funcA(Delta, N, x)
    A = np.concatenate((np.zeros(int((sep-d)/2/Delta)),
                        np.ones(int(d/Delta)), np.zeros(int((sep-d)/2/Delta))))
    A = (np.tile(A, slits))
    A = np.concatenate((np.zeros(int((N-len(A))/2)), A,
                        np.zeros(N-len(A)-int((N-len(A))/2))))
    return(A*fA)

# computes the intensity on the screen a distance D away from an aperture A
# takes the Fourier transform of A using N samples spaced evenly in L and scales
# output accordingly


def FFT(D, N, L, A, Fres=False):
    Delta = (L/N)
    x = np.arange(N)*Delta
    if Fres:
        # quadratic phase term which becomes significant in the Fresnel regime
        A = A*np.exp(1j*(x-N/2*Delta)**2*np.pi/D/lmda)
    phi = fft.fft(A)/N  # normalising wave function
    j = fft.fftshift(fft.fftfreq(N, Delta))
    powerFFT = np.abs(fft.fftshift(phi))**2
    return(D*lmda*j, powerFFT)


L = 5000
lmda = 0.5

#%%
"""
=================================
§1. Far field (Fraunhofer regime)
=================================

----------------
§1.1 Single slit
----------------
"""
# Single slit case
d = 100  # slit width
D = 10**6  # distance from screen
N = 2**11  # number of FFT bins
[y, powerFFT] = FFT(D, N, L, aper(d, N, L))
N = 2**15  # number of FFT bins
[y1, powerFFT1] = FFT(D, N, L, aper(d, N, L))
# standard result for single slit in Fraunhofer regime


def powerTheory(y): return (d/L)**2*np.sinc(d*y/D/lmda)**2


pwerTheory = powerTheory(y)
plt.figure()
plt.xlim(-3, 3)
plt.ylim(0, 1.05)
plt.plot(y/10000, powerFFT/pwerTheory.max(), '.',
         markersize=2, label=r'FFT intensity, $N=2^{11}$')
plt.plot(y1/10000, powerFFT1/pwerTheory.max(), '.',
         markersize=2, label=r'FFT intensity, $N=2^{15}$')
plt.plot(y/10000, pwerTheory/pwerTheory.max(),
         lw=0.7, label='Theoretical intensity')
plt.xlabel(r'$y$/cm')
plt.ylabel(r'Relative Intensity $|\psi|^{2}/|\psi_0|^{2}$')
plt.legend()
plt.savefig('singleSlit.pdf', format="pdf")

#%%
plt.figure()
plt.xlim(-3, 3)
plt.ylim(-0.3, 1.1)

# find and plot errors in the FFT result compared to theoretical result
plt.plot(y/10000, (pwerTheory-powerFFT)/(pwerTheory-powerFFT).max(),
         '-', label=r'FFT intensity, $N=2^{11}$')
plt.plot(y1/10000, (powerTheory(y1)-powerFFT1)/(pwerTheory -
                                                powerFFT).max(), '-', label=r'FFT intensity, $N=2^{15}$')
plt.xlabel(r'$y$/cm')
plt.ylabel(r'Relative Error in $|\psi|^{2}$')
plt.legend()
plt.savefig('singleSlitDiff.pdf', format="pdf")
#%%

"""
-----------------------------
§1.2 Sinusoidal phase grating
-----------------------------
"""
# Finite sinusoidal phase grating
d = 2000
D = 10**7
N = 2**15  # minimum required tested to be 2^11 for accurate result, using 2^15 due to low marginal cost in computational time
m = 8


def funcA(Delta, N, x): return np.exp(
    1j*(m/2)*np.sin(2*np.pi*(x-N/2*Delta)/100))  # s and m hardcoded in
# powerTheory=d/D/lmda*np.sinc(d*y/D/lmda)**2
# lmda*D/100 #expected separation of peaks


[y, powerFFT] = FFT(D, N, L, aper(d, N, L, funcA=funcA))

plt.figure()
plt.xlim(-50, 50)
plt.ylim(0, 1.05)
plt.plot(y/10000, powerFFT/powerFFT.max(), '-', lw=1)
# plt.plot(y/10000,powerTheory*10**5)
plt.xlabel(r'$y$/cm')
plt.ylabel(r'Relative Intensity $|\psi|^{2}/|\psi_0|^{2}$')
plt.savefig('sinSlit.pdf', format="pdf")
#%%
"""
===========================
§1.2.1 Effects of varying s
===========================
import matplotlib.animation as animation

#Finite sinusoidal phase grating
d=2000
D=10**7
N=2**10
s=100
marr=np.linspace(2,25,200)
def init():
    line1.set_data([], [])
    Dtext.set_text('')
    return(line1,Dtext)

def animate(i):
    m=marr[i]
    [y,powerFFT]=FFT(D,N,L,aper(d,N,L,funcA=lambda Delta,N,x: np.exp(1j*(m/2)*np.sin(2*np.pi*(x-N/2*Delta)/100))))
    line1.set_data(y/10000,powerFFT*10**2)
    Dtext.set_text('m=%i' % (marr[i]))
    return (line1,Dtext)

fig = plt.figure()
ax = plt.axes(ylim=(0,3), xlim=(-50,50) ,autoscale_on=False)
ax.set_xlabel(r'$y$/cm')
ax.set_ylabel(r'Relative Intensity $|\psi|^{2}/|\psi_0|^{2}$')
line1, = ax.plot([], [], '-', lw=1)
Dtext = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ani = animation.FuncAnimation(fig, animate, frames=np.arange(1, len(marr)), interval=30, init_func=init)
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='tyl35'), bitrate=1800)
ani.save('m.mp4', writer=writer)
"""
#%%
"""
----------------
§1.3 Double slit
----------------
"""
d = 20  # slit width
D = 10**6  # distance from screen
N = 2**15
sep = 80
slits = 2  # double slit experiment
[y, powerFFT] = FFT(D, N, L, aper(d, N, L, slits=slits, sep=sep))

# standard result for a double slit in Fraunhofer regime
powerTheory = (slits*d/L*np.cos(np.pi*sep*y/D/lmda))**2
powerTheory2 = (2*d/L*np.sinc(d*y/D/lmda))**2

plt.figure()
plt.xlim(-5, 5)
plt.ylim(0, 1.05)
plt.plot(y/10000, powerTheory/powerTheory.max(), '--',
         lw=0.5, label=r'2 slits, width $\delta$')
plt.plot(y/10000, powerFFT/powerTheory.max(), '-', label=r'2 slits, width $d$')
plt.plot(y/10000, powerTheory2/powerTheory.max(),
         '--', lw=0.5, label=r'1 slit, width $d$')
plt.xlabel(r'$y$/cm')
plt.ylabel(r'Relative Intensity, $|\psi|^2/|\psi_0|^{2}$')
plt.legend(loc=1)
plt.savefig('doubleSlit.pdf', format="pdf")
#%%

"""
==============================
§2 Near field (Fresnel regime)
==============================

----------------
§2.1 Single slit
----------------
"""
from scipy import integrate


def c1(u): return ((lambda f1: integrate.quad(f1, 0, u)[0])(lambda x: np.cos(
    np.pi*x**2/2)))  # evaluates C(u) using a scipy.integrate method


c = np.vectorize(c1)


def s1(u): return ((lambda f2: integrate.quad(f2, 0, u)[0])(lambda x: np.sin(
    np.pi*x**2/2)))  # evaluates S(u) using a scipy.integrate method


s = np.vectorize(s1)


def psi(x0, x1, D, lmda): return (c(x1*np.sqrt(2/lmda/D))-c(x0*np.sqrt(2/lmda/D))
                                  ) + 1j*(s(x1*np.sqrt(2/lmda/D))-s(x0*np.sqrt(2/lmda/D)))  # evaluates Ψ


d = 100
D = 5*10**3
N = 2**11

[y, powerFFT] = FFT(D, N, L, aper(d, N, L), Fres=True)
N = 2**15
[y1, powerFFT1] = FFT(D, N, L, aper(d, N, L), Fres=True)
yCornu = np.linspace(-250, 250, 1000)
powerCornu = np.abs(psi(yCornu-d/2, yCornu+d/2, D, lmda))**2*((d/L)**2)/8
# factor of d/L**2 normalises wavefunction incident on the whole aperture, #factor of 8 normalises Cornu spiral intensity from -∞,∞

plt.figure()
plt.xlim(-0.25, 0.25)
plt.ylim(0, 1.1)
plt.plot(y/1000, powerFFT/powerCornu.max(), '.',
         lw=1, ms=1, label=r'FFT, $N=2^{11}$')
plt.plot(y1/1000, powerFFT1/powerCornu.max(), '.',
         lw=1, ms=1, label=r'FFT, $N=2^{15}$')
plt.plot(yCornu/1000, powerCornu/powerCornu.max(),
         '-', lw=0.7, label='Fresnel integral')
plt.xlabel(r'$y$/mm')
plt.ylabel(r'Relative Intensity $|\psi|^{2}/|\psi_0|^{2}$')
plt.legend()
plt.savefig('singleSlitFres.pdf', format="pdf")

(powerCornu.max()-powerFFT1.max())/powerCornu.max()
#%%

"""
-----------------------------
§2.2 Sinusoidal phase grating
-----------------------------
"""
d = 2000
s = 100
D = 5*10**5
N = 2**18  # minimum for sensible result 2^16, using 2^18 due to low marginal cost in computational time

[y, powerFFT] = FFT(D, N, L, aper(d, N, L, funcA=funcA), Fres=True)

plt.figure()
plt.xlim(-25,25)
plt.ylim(0, 1.05)
plt.plot(y/1000, powerFFT/powerFFT.max(), '-', lw=1)
plt.xlabel(r'$y$/mm')
plt.ylabel(r'Relative Intensity $|\psi|^{2}/10^{-3}$')
plt.savefig('sinSlitFres.pdf', format="pdf")

#%%

"""
========================================
§3 Transition between near and far field
========================================

import matplotlib.animation as animation

d=200
D=2**np.linspace(8,17,300)
N=2**15

def init():
    line1.set_data([], [])
    Dtext.set_text('')
    return(line1,Dtext)

def animate(i):
    [y,powerFFT]=FFT(D[i],N,L,aper(d,N,L),Fres=True)
    line1.set_data(y,powerFFT/np.max(powerFFT))
    Dtext.set_text('D=%i' % (D[i]))
    return (line1,Dtext)

fig = plt.figure()
ax = plt.axes(ylim=(0,1), xlim=(-300,300) ,autoscale_on=False)
line1, = ax.plot([], [], '-', lw=1)
Dtext = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ani = animation.FuncAnimation(fig, animate, frames=np.arange(1, len(D)), interval=30, init_func=init)
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='tyl35'), bitrate=1800)
ani.save('transition.mp4', writer=writer)
"""
plt.show()
