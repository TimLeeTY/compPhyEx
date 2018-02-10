# Physics Part II: Computational Physics Exercises

This project contains my solutions for the 3 exercises in the computational physics course, these were all attempted in ```python3``` with a wide usage of the ```numpy``` and ```scipy``` packages, a standard  adopted by the scientific computing community. 

The first exercise was on numerical integration, exploring the use of a standard Gaussian method and a Monte-Carlo method. The implementation of the Gaussian method was straightforward using ```scipy.integrate.quad```, and was used to solve the *Fresnel* integrals. However, this technique scales poorly for higher dimensions which is why a Monte-Carlo method was used for the case of an 8-dimensional integral.

The ```python``` code for this exercise is separated into ```numInt.py``` and ```fresnel.py```. Relavent figures are plotted out in pdf format using ```matplotlib.pyplot```.
