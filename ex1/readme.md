---
title: |
  | **Computational Physics**
  | Exercise 1: Integration and Random Numbers
header-includes:
  - \usepackage{fullpage}
  - \usepackage{xfrac}
  - \usepackage{float}
  - \usepackage{caption}
  - \usepackage{subfig}
  - \newcommand*\diff{\mathop{}\!\mathrm{d}}
author: Tim Lee, tyl35
numbersections: false
---

# Task 1

The code used for this task is contained in §1 of ` numInt.py ` which includes a simple implementation of a Monte-Carlo integration applied to the 8-dimensional integral:

$$\int_0^s \sin (x_0+x_1+\dots +x_7)\diff^8 \pmb{x}. $$

The integrand was evaluated at $N$ samples points within the integral volume with $N\in[2^5,2^6,\dots,2^{21}]$.

\begin{figure}[H]
\captionsetup{width=0.8\textwidth}
\centering
{\includegraphics[width=4.5in]{est.pdf}}
\caption{Estimates of the integral for a range of $N$ plotted on a \texttt{semilogx} scale with errors that correspond to the standard deviation of 25 samples. The analytic answer is also plotted (green)}
\end{figure}

Our estimate of the integral clearly converges to the true value; with up to $N=2^{21}$, our best estimate for the integral was $537.19 \pm 0.02$  (averaged over 25 trials). For a single trial, an answer was obtained using $N=2^{26}$ within 10 seconds which gave the value $537.191 \pm 0.005$. Here the usage of the error estimate is justified as it is shown to be very close to the true error as discussed below.

## Problems

One problem that was encountered when developing the code was a memory limit imposed by the MCS operating system. In particular, the number of elements in a `float64` array could not exceed $\sim 5.5 \times 10^6$, meaning that instead of starting with an array of $8 \times N \times n$ random numbers, a `for` loop over the $n$ trials was implemented to reach larger values of $N$. Fortunately, this proved not to hinder performance given the small values of $n$ used.

# Supplementary task 1

2 methods for evaluating the error were compared; the first involved using equation (3) from the handout:

$$\sigma \approx V\left(\frac{\langle f^2 \rangle-\langle f \rangle ^2}{N}\right)^{\sfrac{1}{2}},$$

which gives an estimate of the error, and the second method involved performing the Monte-Carlo integration 25 times and taking the standard deviation of results for each $N$ which is justified as we expect the random numbers sampled to be uncorelated. The code for this section is included in §2 of `numint.py`.

\begin{figure}[H]
\captionsetup{width=0.8\textwidth}
\centering
{\includegraphics[width=4.5in]{err.pdf}}
\caption{Errors of the integral, $\Delta$, for a range of $N$ found using (I) the standard deviation of the mean (blue) and (II) $\sigma$ from equation (3) (orange), along with a `line of best fit' for each, plotted on a \texttt{loglog} scale.}
\end{figure}

It is clear that the errors derived from either method agree to a large extent suggesting equation (3) is a good estimate of the error. To verify that the errors scaled as $N^{-\sfrac{1}{2}}$, we determined the gradient of the $\ln\Delta$ against $\ln N$ and found that they were within $-0.50\pm0.02$ as expected.

# Task 2

The primary goal of this task is to tackle the Fresnel integral:

$$\int_0^u \cos\left(\frac{\pi x^2}{2}\right)+i\sin\left(\frac{\pi x^2}{2}\right)\diff x.$$

The code used to tackle this is included in §1 of `fresnel.py`, the real part of the integral is evaluated by `c(u)` and the imaginary part by `s(u)`.

\begin{figure}[H]
\captionsetup{width=0.8\textwidth}
\centering
{\includegraphics[width=4.5in]{FresnelSpiral.pdf}}
\caption{A plot of the Cornu spiral on the complex plane.}
\end{figure}

# Supplementary task 2

By scaling the x coordinate accordingly, `c(u)` and `s(u)` are used to evaluate the relative amplitude and phase of the wavefunction at the screen. The code used to tackle this is included in §2 of `fresnel.py`. 

\begin{figure} \captionsetup{width=0.8\textwidth}
\centering
{\includegraphics[width=4.5in]{FresnelSlitAmp.pdf}}
\caption{The relative intensity of the diffraction pattern on the screen along $x$ for the three different slit widths $D$.}
\end{figure}

\begin{figure}
\captionsetup{width=0.8\textwidth}
\centering
{\includegraphics[width=4.5in]{FresnelSlitPha.pdf}}
\caption{The relative phase of the diffraction pattern on the screen along $x$ for the three different slit widths $D$.}
\end{figure}
