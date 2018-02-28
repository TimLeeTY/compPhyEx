---
title: |
  | **Computational Physics**
  | Exercise 3A: Diffraction by the FFT
header-includes:
  - \usepackage{fullpage}
  - \usepackage{xfrac}
  - \usepackage{float}
  - \usepackage{subfig}
  - \usepackage{siunitx}
  - \usepackage{caption}
  - \newcommand*\diff{\mathop{}\!\mathrm{d}}
  - \DeclareMathOperator{\sinc}{sinc}
author: Tim Lee, tyl35
numbersections: false
---

# Task 1

Our aim is to find the wavefunction of the diffraction pattern, $\psi(y)$, as a function of $y$, by integrating over the aperture. This is solved by appropriately applying the Fourier transform to the aperture function $A(x)$ which can be can be written in the discretised form:

$$\psi(y)\propto \Delta \sum_{m=0}^{N-1} \underbrace{\exp\left(\frac{ikx_m^2}{2D}\right)A(x_m)}_{A'(x_m)} \exp\left(\frac{-ikyx_m}{D}\right),$$ {#eq:1}

where we define $x_m=(m-\sfrac{N}{2})\Delta$. In the Fraunhofer limit, $A'(x)\to A(x)$. Comparing this to the definition of a general discrete Fourier transform:

$$H_j=\sum_{m=0}^{N-1}h_m\exp\left(\frac{2\pi i \,m j}{N}\right)$$

and defining the Fourier frequency: $f_j=j/(N\Delta)$, the discrete Fourier transform of the aperture function $A(x_m)$ can be rewritten as:

$$\psi(y_j)=\sum_{m=0}^{N-1}A(x_m)\exp\left(2\pi i \,m\Delta f_j\right).$$ {#eq:2}

Comparing @eq:1 with @eq:2 tells us that we can map $f_j$, the Fourier frequency, to $y$, the desired displacement by $y_j=f_jD\lambda$.

The code used to tackle this exercise is all included in `FFT.py`. 2 main functions were written, `aper` and `FFT`, the former takes arguments that describe the aperture (like slit width and number of slits) and returns the aperture function as an array, `A` as a function of $x$ with increments of $\Delta$. `FFT` applies the fast Fourier transform (included in the `numpy.fft` package) and returns the power spectrum as an array along with appropriately scaled $y$ values.

\begin{figure} [H]
  \centering
	\captionsetup{width=0.95\textwidth}
  \subfloat[Relative intensity]{
    {\includegraphics[width=3in]{singleSlit.pdf}}}\quad
  \subfloat[Error in intensity]{
    \includegraphics[width=3in]{singleSlitDiff.pdf} \label{singleError}}
  \caption{Plot of the intensity pattern of a single slit aperture as a function of distance on the screen (a), obtained from a FFT, and their accompanying errors (b) compared to the theoretical intensity. 2 different numbers of samples were compared: $N=2^{11}$ (blue) and $N=2^{15}$ (orange).}
  \label{singleSlit}
\end{figure}

The power spectrum obtained for a single slit patter was then computed and plotted (figure \ref{singleSlit}). The output could then be compared to an analytic solution of the full Fourier transform:

$$|\psi_\mathrm{single}|^2=\left(\frac{d}{L}\right)^2 \sinc^2\left(\frac{d}{D\lambda}y\right),$$ {#eq:singleFunc}

where the factor $d/L$ is included to ensure the incident wavefunction on the whole aperture is normalised.

Figure \ref{singleError} shows that the result obtained from the FFT agrees with the theoretical solution for the most part, with a notable exception around the region of $y=0$ where the FFT appears to fall slightly short of maximum intensity. This appears to be a problem with undersampling as the error greatly diminishes with more samples (the errors are $\sim 50$ times smaller with a $2^4$ times increase in $N$ from $2^{11}$ to $2^{15}$).

## Further work: double slits

The aperture function `aper` was modified to be able to produce multiple slit apertures with some separation. While still far from a completely general aperture function, this allowed us to test whether the FFT produced reasonable patterns for the well known multiple slits cases.

In particular, the double slits case was chosen due to its simple analytic solution. The pattern for 2 slits with width $d$ separated by $a$ was plotted along with the analytic solutions of:

1. single slit with width $d$ (@eq:singleFunc);
2. two infinitely narrow slits separated by $a$:
  $$|\psi_\mathrm{double}|^2\propto\cos^2\left(\frac{\pi a }{D\lambda}y\right).$$

The intensity patterns are plotted in §1.3 of `FFT.py` and are shown in figure \ref{doubleSlit}.

\begin{figure}[H]
\captionsetup{width=0.9\textwidth}
\centering
{\includegraphics[width=3.5in]{doubleSlit.pdf}}
\caption{Plot of the intensity pattern of finite double slits (orange) as a function of distance on the screen, $y$ obtained using $N=2^{15}$. This is compared to the analytic solutions of a single slit (green dashed) and double slits (blue dashed).}
\label{doubleSlit}
\end{figure}

# Task 2

The code used for this task is included in §1.2 of `FFT.py`. The single slit from the previous section is modified into a sinusoidal phase grating by multiplying the aperture function, $A(x)$ by $e^{i\phi(x)}$ where:

$$\phi(x)=\frac{m}{2}\sin\left(\frac{2\pi x}{s}\right)$$

through the `funcA` parameter of the `FFT` function. The resultant intensity is shown in figure \ref{sinSlit}.

\begin{figure}[H]
\captionsetup{width=0.9\textwidth}
\centering
{\includegraphics[width=3.5in]{sinSlit.pdf}}
\caption{Plot of the intensity pattern of a sinusoidal phase grating as a function of distance on the screen, $y$ obtained using $N=2^{15}$.}
\label{sinSlit}
\end{figure}

As the aperture is formed by a multiplication of an infinite sinusoidal grating with a finite slit, we expect the pattern to be a convolution of the single slit pattern with that of an infinite sinusoidal grating. This is confirmed by the fact that each of the spikes (spaced at regular $\SI{5}{\centi\meter}$ intervals) appear to form the single slit pattern upon closer inspection. To better understand this we turn to the aperture integral of an infinite sinusoidal phase grating:

$$\Psi(y)\propto\int\diff x \, \exp\left(\frac{im}{2}\sin\left(\frac{2\pi x}{s}\right)\right) \exp\left(\frac{-ikyx}{D}\right).$$

While the integrand can not be readily evaluated without the use of Bessel functions, it can be appreciated that the integral will only be non-zero when the wavelength of the 2 terms match (by the orthogonality of sine and cosine.), giving rise to maxima at regular intervals. Further treatment of the integral reveals that the amplitude of the $n^{\mathrm{th}}$ peak from the centre is $\propto J_n(m/2)$ where $J_n$ is the $n^{\mathrm{th}}$ order Bessel function of the first kind, while their separation is given by $\lambda D/s$.

To further investigate this, an animation was set up using the `matplotlib.animation` package in §1.2.1 of `FFT.py`, which shows how the form of the intensity pattern changes with different $m$. As expected, it showed that varying $m$ had no effect on the separation of peaks, only their relative intensities. The result was outputted as an `.mp4` file and was therefore not submitted as part of this exercise.

# Supplementary Task

The code used for this task is included in §2 of `FFT.py`. The implementation of Fresnel effects was done simply by multiplying the aperture function by the quadratic exponential $\exp\left(ikx^2/2D\right)$ (see @eq:1), returning to the non-approximated form of the aperture integral.

\begin{figure}
\captionsetup{width=0.9\textwidth}
\centering
{\includegraphics[width=3.5in]{singleSlitFres.pdf}}
\caption{Plot of the intensity pattern of a single slit aperture as a function of distance on the screen in the Fresnel regime, obtained using a FFT for $N=2^{11}$ (blue) and $N=2^{15}$ (orange) compared to the theoretical intensity (green).}
\label{singleFres}
\end{figure}

For the single slit case (shown in figure \ref{singleFres}), we compared the result from the FFT to the theoretical solution from Exercise 1. The result clearly shows that the accuracy of the FFT improves for larger $N$. Furthermore, the FFT with smaller $N$ appears to break symmetry about $x=0$, assumed to be an artifact of the dicretisation of the aperture function. Perhaps surprisingly, compared to the Fraunhofer case, no increase in sampling rate was required to achieve a similar level of accuracy for a single slit (both had a ~0.1% error with $N=2^{15}$).

\begin{figure}[H]
\captionsetup{width=0.9\textwidth}
\centering
{\includegraphics[width=3.5in]{sinSlitFres.pdf}}
\caption{Plot of the intensity pattern of a sinusoidal phase grating as a function of distance on the screen, $y$, in the Fresnel regime obtained using $N=2^{18}$.}
\label{sinSlitFres}
\end{figure}

The result for the sinusoidal phase grating is perhaps more surprising, while not entirely dissimilar to the far field case (figure \ref{sinSlit}), the patten does appear to be noisier with a distinct periodic pattern convoluted with the original delta peaks of the Fraunhofer case. Perhaps more curiously, it seems to break inversion symmetry about $y=0$ as the peaks on the left are noticeably greater than their counterparts on the right. This can be reconciled with the asymmetry in the phase grating which has an odd parity about $x=0$. 

Unlike the single slit case, the FFT performed in this case required $N=2^{16}$ before the pattern converged to the one shown in figure \ref{sinSlitFres} (no noticeable difference for $N$ up to $2^{22}$), a factor of $2^5$ greater than the minimum of $2^{11}$ for the Fraunhofer case. 

Finally, again as a small extra challenge, a second animation was created to show the intensity patterns in the regime where $D \sim d^2/2$ for a single slit. This was an attempt to reconcile the two extreme limits and was done mainly out of personal curiosity. See §3 of `FFT.py`.
