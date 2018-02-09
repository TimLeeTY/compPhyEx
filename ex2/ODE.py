"""
=========================================
Ex.2 Solving the ODE of a simple pendulum
=========================================
Solve equation: θ''=-(g/l)*sin(θ)-q*θ'+F*sin(Ω*t)
Define: θ=θ, ⍵=θ'
The coupled first order ODES are:
⍵=θ'
⍵'=-(g/l)*sin(θ)-q*⍵+F*sin(Ω*t)
Apply 4th order Runge-Kutta
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


def dy(y,t): return( np.array([y[1], -1*np.sin(y[0])-q*y[1]+F*np.sin(Omg*t)])) #returns the derivatives of θ and ⍵, F and q referenced from global scope (bad practice but neater code)

def RK4step(dy,t,y,dt): #More 'standard' implementation of 4th order Runge-Kuta method
    k1=dt*dy(y,t)
    k2=dt*dy(y+k1/2,t+dt/2)
    k3=dt*dy(y+k2/2,t+dt/2)
    k4=dt*dy(y+k3  ,t+dt)
    return(k1/6+k2/3+k3/3+k4/6)

def RK4(dy): #Implementation of 4th order Runge-Kuta method using Lambda notation
    return( lambda t,y,dt:( #takes arguments of t (time), y (the current [θ,⍵]) and dt (the timestep)
            lambda k1:(
            lambda k2:(
            lambda k3:(
            lambda k4: k1/6+k2/3+k3/3+k4/6 #returns the vector [Δθ(t),Δ⍵(t)]
            )(dt*dy(y+k3  ,t+dt))
            )(dt*dy(y+k2/2,t+dt/2))
            )(dt*dy(y+k1/2,t+dt/2))
            )(dt*dy(y,t))
            )

def solveRK(RK4,dy,dt,Tf,y0,y1=0): #performs integration given dy (the derivatives), dt (the timestep), Tf (the final time) and y0,y1 (the initial conditions)
    t,y=0.,np.array([y0,y1])         #initial conditions
    finalRK=np.zeros((int(Tf/dt),2)) #initialises array, final holds results as a list of  [θ(t),⍵(t)] values
    for i in range(len(finalRK)):
        finalRK[i]= y
        #t,y=t+dt,y+RK4step(dy,t,y,dt)
        t,y=t+dt,y+RK4(dy)(t,y,dt)
    return(finalRK) #returns an 2D array with rows of [θ(t),⍵(t)] at each timestep

#%%
q,Omg,F=0,2/3,0               #setting parameters
dt=0.05

"""
------------------------
§1. Period vs. amplitude
------------------------
Investigating how the period, T, is affected by θ_0
"""

def period(tarr,y): #finds the avearge time between when y changes sign
    return( lambda sgn: (
            lambda T: np.array([np.mean(T[1:]-T[:-1]),np.std(T[1:]-T[:-1])]) #returns time period and associated error
    )((tarr[:-1])[(sgn[1:]-sgn[:-1])!=0]*2) #time stamps for when y changes sign (*2 for actual period)
    )(sgn=np.sign(y)) #sign of y

Tf=600 #average time period over ~100 oscillations
tarr=np.arange(0.0, Tf, dt)
y0=np.linspace(0.01,np.pi-0.01,100)
parr=np.array([period(tarr,integrate.odeint(dy,[i,0],tarr)[:,1]) for i in y0]) #more efficient LSODA method used here

#plotting results
fig, ax = plt.subplots()
plt.xlabel(r"$\theta_0$")
plt.ylabel(r"$T$/s")
plt.xlim([0,np.pi])
plt.ylim([5,30])
ax.plot(y0,parr[:,0])
fig.savefig("TvsAmp.pdf",format="pdf")

#%%
"""
-------------------
§2. Varying q and F
-------------------
Investigating how θ and ⍵ vary with different values of q and F
(θ_0 kept at 0.2 throughout)
"""
tD=np.pi*3
Tf=100
tarr=np.arange(0.0, Tf, dt)
F,qarr=0,[5,2,0.5] #keeping F at 0, varying q

#plotting θ in for loop on one figure
fig, ax = plt.subplots()
plt.xlabel(r'$t/T_\mathrm{0}$')
plt.ylabel(r'$\theta (t)$/rad')
plt.xlim([0,6])
for i in qarr:
    q=i
    ax.plot(tarr/2/np.pi,solveRK(RK4,dy,dt,Tf,0.2)[:,0],label='q={:.1f}'.format(q),lw=0.8)
ax.legend()
fig.savefig('q_theta.pdf',format="pdf")

#plotting ⍵ in for loop on one figure
fig, ax = plt.subplots()
plt.xlabel(r'$t/T_\mathrm{0}$')
plt.ylabel(r'$\omega (t)$/rad s$^{-1}$')
plt.xlim([0,6])
for i in qarr:
    q=i
    ax.plot(tarr/2/np.pi,solveRK(RK4,dy,dt,Tf,0.2)[:,1],label='q={:.1f}'.format(q),lw=0.8)
ax.legend()
plt.gcf()
fig.savefig('q_omega.pdf',format="pdf")


q,Farr=0.5,[0.5, 1.2, 1.44, 1.465] #keeping q fixed at 0.5 while varying F

#plotting θ in for loop on one figure
fig, ax = plt.subplots()
plt.xlabel(r'$t/T_\mathrm{D}$')
plt.ylabel(r'$\theta (t)$/rad')
plt.xlim([0,8])
for i in Farr:
    F=i
    ax.plot(tarr/tD,solveRK(RK4,dy,dt,Tf,0.2)[:,0],label='F={:.3f}'.format(F),lw=0.8)
ax.legend()
fig.savefig('F_theta.pdf',format="pdf")

#plotting ⍵ in for loop on one figure
fig, ax = plt.subplots()
plt.xlabel(r'$t/T_\mathrm{D}$')
plt.ylabel(r'$\omega (t)$/rad s$^{-1}$')
plt.xlim([0,8])
plt.ylim([-3,3])
for i in Farr:
    F=i
    ax.plot(tarr/tD,solveRK(RK4,dy,dt,Tf,0.2)[:,1],label='F={:.3f}'.format(F),lw=0.8)
ax.legend(fontsize=8)
fig.savefig('F_omega.pdf',format="pdf")

#%%

"""
-----------
§3. θ vs. ⍵
-----------
Investigating the form of θ vs. ⍵ with different values of F
Shows example of period doubling and quadrupling
"""
Tf=50*tD
tarr=np.arange(0.0, Tf, dt)

q,Farr=0.5,[0.5, 1.2, 1.44, 1.465] #keeping q fixed at 0.5 while varying F
title=dict(zip(Farr,['Small angle','Chaotic','Period doubling','Period quadrupling']))
#plotting different F values on separate figures
for i in Farr:
    F=i
    plt.figure()
    plt.xlabel(r'$\theta (t)$/rad')
    plt.xlim([-np.pi,np.pi])
    plt.ylabel(r'$\omega (t)$/rad s$^{-1}$')
    y=integrate.odeint(dy,[0.2,0],tarr)
    plt.plot((y[:,0]+np.pi)%(2*np.pi)-np.pi,y[:,1],'.',markersize=1)
    #centring θ around 0 and in the range [-π,π]
    plt.title('{}'.format(title[F]))
    plt.gcf()
    plt.savefig('y0vy1_F{:.0f}.pdf'.format(F*1000),format="pdf")
#%%

"""
--------------------
§4. 'Chaotic motion'
--------------------
Investigating the sensitivity of motion to initial conditions for F=1.2
"""
F,q=1.2,0.5
Tf=500
tarr=np.arange(0.0, Tf, dt)
y0=[0.2,0.20001,0.20000001]  #initial θ_0 vary by only 1/20000

#plotting both graphs on the same axis
fig, ax = plt.subplots()
plt.xlim(0,Tf/tD)
plt.xlabel(r'$t/T_\mathrm{D}$')
plt.ylabel(r'$\theta (t)$/rad')
for i in y0:
    ax.plot(tarr/tD,solveRK(RK4,dy,dt,Tf,i)[:,0], label=r"$\theta_0={}$".format(i))
ax.legend()
fig.savefig('chaotic.pdf',format="pdf")


#%%
"""
-------------------
§5. Plotting Energy
-------------------
Comparing 4th order Runge-Kuta method implemented above with the standard scipy
LSODA implementation
Computing time is substantial (the R-K implementation in particular, to save on
time, results are plotted separately by 'energyPlot.py' as to avoid reevaluating
the integral every time I wish to change the plot.
"""
q,Omg,F=0,2/3,0               #setting parameters
dt=0.05                       #step size
Tf=int(np.pi*2*10000)              #final time (10000 natural oscilaltions)
tarr=np.arange(0.0, Tf, dt)
y0=0.01                       #initial displacement
t,y=0.,np.array([y0,0.])      #initial conditions

finalRK=solveRK(RK4,dy,dt,Tf,y0)       #my implementation of 4th order R-K, slower by factor of ~10
finalLSODA=integrate.odeint(dy,y,tarr) #default integrate.odeint implementation of LSODA
trunc=np.arange(len(tarr))%(2**8)==0 #keeping every 256 entries to decrease file size
tarr,finalRK,finalLSODA=tarr[trunc],finalRK[trunc],finalLSODA[trunc]
np.savetxt('energyDiffOut.csv', np.stack((tarr,finalRK.T[0],finalRK.T[1],finalLSODA.T[0],finalLSODA.T[1])).T, delimiter=',')

#%%

"""

===================================
Appendix - Animating chaotic motion
===================================
Commented out as it does not directly relate to the task and provides no
particular insight into the physics.

import matplotlib.animation as animation

F,q=1.2,0.5
Tf=300
tarr=np.arange(0.0, Tf, dt)
finalRK1=solveRK(RK4,dy,dt,Tf,0.2)
finalRK2=solveRK(RK4,dy,dt,Tf,0.20001)
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return(line1,line2)

def animate(i):
    tempx=[0,np.sin(finalRK1[3*i][1])]
    tempy=[0,-np.cos(finalRK1[3*i][1])]
    tempx1=[0,np.sin(finalRK2[3*i][1])]
    tempy1=[0,-np.cos(finalRK2[3*i][1])]
    line1.set_data(tempx,tempy)
    line2.set_data(tempx1,tempy1)
    return (line1,line2)

fig = plt.figure()
ax = plt.axes(xlim=(-1.2,1.2), ylim=(-1.2,1.2),yticks=[],xticks=[] ,autoscale_on=False)
ax.set_aspect('equal')
line1, = ax.plot([], [], 'o-', lw=2)
line2, =ax.plot([], [], 'o-', lw=2)

ani = animation.FuncAnimation(fig, animate, frames=np.arange(1, int(len(finalRK1)/3)), interval=20, init_func=init)
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='tyl35'), bitrate=1800)
ani.save('pendulum.mp4', writer=writer)

plt.show()
"""
