---
title: |
  | **Computational Physics**
  | Exercise 2: Ordinary Differential Equations
header-includes:
  - \usepackage{fullpage}
  - \usepackage{xfrac}
  - \usepackage{float}
  - \usepackage{siunitx}
  - \newcommand*\diff{\mathop{}\!\mathrm{d}}
author: Tim Lee, tyl35
numbersections: false
---

# Task 1

The code used to tackle this task is contained in ```ODE.py``` and  ```energyPlot.py```. 2 methods for solving the coupled ODEs were used, first a 4th order Runge-Kuta method is implemented as the function ```RK4``` which finds the corresponding change in $\theta$ and $\omega$ for each time step. The second method involved using the ```integrate.odeint``` function from the ```scipy``` package which, according to online documentation, uses an implementation of [LSODA]. In terms of computation time, the method built in to ```scipy``` easily out-performed this particular implementation of the Runge-Kuta method by around a factor of 10.

[LSODA]:http://www.oecd-nea.org/tools/abstract/detail/uscd1227

## Energy Conservation

The total energy of the system (kinetic and potential) is calculated and used to test how well either method conserves energy in the undriven and undamped case for up to $\sim 10^4$ oscillations. The code for this is in ```energyPlot.py``` and §5 of ```ODE.py```. Due to long computation times, the results are written into ```energyDiffOut.csv``` which is then read and plotted in ```energyPlot.py```.

\begin{figure}[H]
\captionsetup{width=0.9\textwidth}
\centering
{\includegraphics[width=3.5in]{energyDiff.pdf}}
\caption{The fractional change in energy over time is plotted for the 2 methods on separate vertical axes. The left axis (blue) corresponds to the Runge-Kuta method and the right (orange) corresponds to the LSODA. The time is in units of $T_\mathrm{0}$, the natural period.}
\end{figure}

While the Runge-Kuta method is inefficient time-wise, it produces a more accurate result, outperforming the LSODA method by a factor of $\sim 500$ in terms of energy conservation. This demonstrates the trade-off between computation time and accuracy.

## Period vs. amplitude

The period of oscillations is found by determining the average time between when $\omega$ changes sign. This is handled by passing the angular velocity vector into ```period```. The period is then found for a range of initial displacements, $\theta_0\in [0,\pi]$ for the undamped and undriven case. The code for this is in §1 of ```ODE.py```.

\begin{figure}[H]
\captionsetup{width=0.9\textwidth}
\centering
{\includegraphics[width=3.5in]{TvsAmp.pdf}}
\caption{The average time period of oscillations, $T$, plotted against the initial displacement, $\theta_0$.}
\end{figure}

As expected, in the limit of small displacements, the period of oscillations agrees with the theoretical result using small-angle approximations with $T\approx 2\pi$. For the case where initial displacement, $\theta_0=\pi/2$, it was found that the period was $7.42\pm \SI{0.04}{\second}$.

# Task 2

We now turn towards the effects of adding damping and a driving force to the pendulum.  The code for this is in §2 of ```ODE.py```.

## Effects of damping

For small oscillations, the three regimes of damping are

$$ \frac{2q}{\sqrt{g/l}}=\zeta
\begin{cases}
<1 & \text{underdamped}\\
=1 & \text{critically damped}\\
>1 & \text{overdamped}
\end{cases}$$

As we took $g=l=1$, critical damping occurs when $q=2$. Hence we chose to plot $\theta$ and $\omega$ for $q=0.5,2,5$. This is done in the first half of §2 of ```ODE.py```.

\begin{figure}[H]
  \centering
	\captionsetup{width=0.95\textwidth}
  \subfloat[$\theta$ vs. time]{\includegraphics[width=3in]{q_theta.pdf}}\quad
  \subfloat[$\omega$ vs. time]{\includegraphics[width=3in]{q_omega.pdf}}
	\caption{Plots of $\theta$ and $\omega$ as functions of $t$, for the underdamped (green), critically damped (orange) and overdamped (blue) cases. The time is in units of $T_\mathrm{0}$, the natural period.}
\end{figure}

For the underdamped case ($q=0.5$), oscillations are still observed, but the amplitude is greatly attenuated over 3 cycles after which it remains at the equilibrium position. For the critically damped case ($q=2$), the displacement of the pendulum returns to equilibrium position very rapidly without crossing $\theta=0$. Any further increase in damping only lengthens the time required for the pendulum to return to equilibrium as demonstrated by the overdamped case ($q=5$). This can be explained by the form of the angular velocities which shows a larger peak in the critically damped case than the overdamped case, suggesting a more rapid decay.

## Effects of a driving force

Next, we fix $q=0.5$ and vary the amplitude of the driving force (with frequency $\Omega_\mathrm{D}= \SI{2/3}{\radian \second^{-1}}$) with $F=0.5, 1.2, 1.44, 1.465$. This is done in the second half of §2 of ```ODE.py```.

\begin{figure}[H]
  \centering
	\captionsetup{width=0.95\textwidth}
  \subfloat[$\theta$ vs. time]{
    \includegraphics[width=3in]{f_theta.pdf}
    \label{f_theta}
  }\quad
  \subfloat[$\omega$ vs. time]{
    \includegraphics[width=3in]{f_omega.pdf}
    \label{f_omega}
  }
	\caption{Plots of $\theta$ and $\omega$ as functions of $t$, for three different amplitudes of driving force. The time is in units of $T_\mathrm{D}$, the period of the driving force.}
  \vspace{10px}
  \label{varyF}
\end{figure}

With a small driving force ($F$ up to $\sim 0.7$), the pendulum is reasonably well behaved and displays a sinusoidal displacement. This is exemplified in figure \ref{varyF} with the $F=0.5$ case (blue) which shows that both $\theta$ and $\omega$ follow a sinusoidal pattern. As we would expect, the period is close to the driving period of $\SI{3 \pi}{\second}$, and using the previously described method, the period was found to be $9.41\pm\SI{0.18}{\second}$ (for $\sim 10^3$ oscillations).

For larger values of $F$, the pendulum begins to behave more erratically with displacements that deviate much further from equilibrium. For the particular case where $F=1.2$ in figure \ref{varyF}, neither graphs of $\theta$ nor $\omega$ show signs of regular behavior. Hence there is no meaningful way to assign a period for $F=1.2$. This is more evident in the next section with plots of $\omega$ against $\theta$.

For even larger values of $F$, we find 2 special cases where $F=1.44, 1.465$ where the period doubles and quadruples respectively. While seemingly identical in figure \ref{f_theta}, the subtle differences are more apparent in figure \ref{f_omega} from which we can see that while the $F=1.44$ graph repeats itself every 2 periods, slight variations for the $F=1.465$ case means it only repeats itself after 4 periods. This is not immediately apparent in these graphs but is made more evident in the following sections with graphs of $\omega$ against $\theta$.

## $\theta$ vs. $\omega$ graphs

To better demonstrate the changes in period described in the previous section, we plotted graphs of $\omega$ against $\theta$. The code for this is in §3 of ```ODE.py```. The graphs shown in figure \ref{omegaVStheta} are for durations of $50\,T_\mathrm{D}$.


\begin{figure}
  \centering
	\captionsetup{width=0.95\textwidth}
  \subfloat[$F=0.5$]{
    \includegraphics[width=3in]{y0vy1_F500.pdf}
    \label{f=0.5}
  }\quad
  \subfloat[$F=1.2$]{
    \includegraphics[width=3in]{y0vy1_F1200.pdf}
    \label{f=1.2}
  }\\
  \subfloat[$F=1.44$]{
    \includegraphics[width=3in]{y0vy1_F1440.pdf}
    \label{f=1.44}
  }\quad
  \subfloat[$F=1.465$]{
    \includegraphics[width=3in]{y0vy1_F1465.pdf}
    \label{f=1.465}
  }
	\caption{Plots of $\omega$ against $\theta$ over $50 T_\mathrm{D}$ for 4 different amplitudes of driving force.}
  \label{omegaVStheta}
\end{figure}

Figure \ref{f=0.5} shows that, for small driving amplitudes, the system follows a very regular pattern which repeats itself every $T_\mathrm{D}$ with one tidy loop (ignoring initial transient effects). When $F=1.2$ (figure \ref{f=1.2}), the system shows no signs of any regular pattern and seemingly never repeats itself over $50\,T_\mathrm{D}$, further reinforcing our previous assertion that no meaningful period can be assigned to it.

For the 2 special cases, we can clearly see that while figure \ref{f=1.44} has 2 distinct loops, figure \ref{f=1.44} shows 4 (2 outer and 2 inner), less distinct loops. Again, this confirms our previous statements about the period doubling for $F=1.44$ and quadrupling for $F=1.465$.

# Supplementary task 1

Finally we investigate the sensitivity of resulting motion to initial case in the 'erratic' case where $F=1.2$. The code for this is in §4 of ```ODE.py```. $\theta$ is plotted as functions of time for 3 different initial displacements $\theta_0=\SI{0.2}{\radian}$, $\theta_0=\SI{0.20001}{\radian}$ and $\theta_0=\SI{0.20000001}{\radian}$, a difference of 1 in $2\times 10^3$ and 1 in $2\times 10^6$ respectively,

\begin{figure}[H]
\captionsetup{width=0.9\textwidth}
\centering
{\includegraphics[width=3.5in]{chaotic.pdf}}
\caption{Displacement as functions of time for initial conditions $\theta_0=\SI{0.2}{\radian}$ (blue), $\theta_0=\SI{0.20001}{\radian}$ (orange) and $\theta_0=\SI{0.20000001}{\radian}$ (green). The time is in units of $T_\mathrm{D}$, the period of the driving force.}
\end{figure}

While the solutions initially seem identical, at around $10 \, T_\mathrm{D}$, they diverge dramatically and become increasingly dissimilar. In this case, we have used the more accurate Runge-Kuta method which has shown to be correct to within 1 in $10^8$ over $10000$ natural oscillations, suggesting that this is not a fault in the integration and that the system is indeed very sensitive to initial conditions.

As we would expect, smaller differences in initial conditions meant it took longer for the displacement to deviate. However, this deviation time seems to increase very slowly even for significant changes in $\Delta\theta_0$. This has been tested for differences in $\theta_0$ down to 1 in $2\times 10^{10}$ after which point errors in the integration method is expected to dominate.

Finally, as an extra challenge, the pendulums with $\theta_0=\SI{0.2}{\radian}$ and $\theta_0=\SI{0.20001}{\radian}$ were animated in the final section of the code which utilises the ```matplotlib.animation``` package to show the pendulum as it evolves over time and is outputted as an ```.mp4``` file. However, it has proven to be computationally intensive and has not provided any meaningful physical insight and so has been left commented out.
