# 1D_Advection-Flux_Reconstruction_Method

Solves the 1D linear and non-linear advection equation using the Flux Reconstruction Discontinuous Galerkin Method of Huynh(2007).

Enables computations using polynomials upto 9th degree.
Time-stepping using three stage Strongly Stability Preserving Runge Kutta method(RK3SSP).
Roe's Riemann solver is employed to compute interface flux.

Periodic boundary condition.
Initial conditions include linear case with Gaussian profile, u(x,0)=exp(-20x^2) 
and non-linear case with sinusoidal profile, u(x,0)=sin (pi x).
