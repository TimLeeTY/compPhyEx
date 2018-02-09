"""
=================================================
Ex.1 Task 2: Fresnel integral and the single slit
=================================================
Fresnel integral:
∫cos(π*x^2/2) + i*sin(π*x^2/2) dx
Integrate real and imaginary parts separately
"""
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def c1(u): return ((lambda f1: integrate.quad(f1, 0, u)[0]) (lambda x: np.cos(np.pi*x**2/2))) #evaluates C(u) using a scipy.integrate method
c=np.vectorize(c1)
def s1(u): return ((lambda f2: integrate.quad(f2, 0, u)[0]) (lambda x: np.sin(np.pi*x**2/2))) #evaluates S(u) using a scipy.integrate method
s=np.vectorize(s1)

"""
----------------
§1. Cornu Spiral
----------------
Plots the Fresnel spiral on the complex plane using c1 and s1 from above
"""

x = np.linspace(-6,6,1000)
fig, ax = plt.subplots()
ax.plot(c(x),s(x),lw=0.8)
ax.set_aspect('equal') #equalises the x and y axis to make the spiral circular
ax.set_xlim([-0.8,0.8])
ax.set_xticks(np.arange(-0.8,0.9,0.2))
ax.set_ylim([-0.8,0.8])
ax.set_xlabel(r'$C(u)$')
ax.set_ylabel(r'$S(u)$')
fig.savefig("FresnelSpiral.pdf",format="pdf")
#%%
"""
---------------------------
§2. Single slit diffraction
---------------------------
Scales the Fresnel integrals and finds the amplitude and phase of light on a
screen formed from the diffraction of a single slit.
"""

xl=20 #extent of pattern on screen
x=np.linspace(-xl,xl,200)
d,D,lmda=10,np.array([30,50,100]),1 #set parameters of diffraction
psi= lambda x0,x1,D,lmda: (c(x1*np.sqrt(2/lmda/D))-c(x0*np.sqrt(2/lmda/D)))+ 1j*(s(x1*np.sqrt(2/lmda/D))-s(x0*np.sqrt(2/lmda/D))) #evaluates Ψ

phiarr=[psi(x-d/2,x+d/2,i,lmda) for i in D] #find complex amplitude of light for various values of D
fig, ax = plt.subplots()
for i in range(len(phiarr)):
    ax.plot(x,np.abs(phiarr[i]), label="D=%i cm"%D[i],lw=0.8) #plot relative amplitude of light
ax.set_xlim([-xl,xl])
ax.set_ylim([0,2])
plt.xlabel(r'$x\,$/cm')
plt.ylabel(r'Amplitude, $|\psi|_{\mathrm{rel}}$')
ax.legend()
fig.savefig("FresnelSlitAmp.pdf",format="pdf")
#%%

xl=9
fig, ax = plt.subplots()
for i in range(len(phiarr)):
    ax.plot(x,np.angle(phiarr[i]), label="D=%i cm"%D[i],lw=0.8) #plot relative phase of light
ax.set_xlim([-xl,xl])
ax.set_ylim([0,3])
plt.xlabel(r'$x\,$/cm')
plt.ylabel(r'Phase, $\theta\,$/rad')
ax.legend(loc=9)
fig.savefig("FresnelSlitPha.pdf",format="pdf")
plt.show()
