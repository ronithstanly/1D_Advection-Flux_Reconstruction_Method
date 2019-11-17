# 1D_Advection_using_Flux_Reconstruction_Method

Solves the 1D linear and non-linear advection equation using the Flux Reconstruction (Discontinuous Galerkin Spectral Element) Method of Huynh(2007).

Enables computations using polynomials upto 9th degree.
Time-stepping using three stage Strongly Stability Preserving Runge Kutta method(RK3SSP).
Roe's Riemann solver is employed to compute interface flux.

Periodic boundary condition.
Initial conditions include linear case with Gaussian profile, u(x,0)=exp(-20x^2) 
and non-linear case with sinusoidal profile, u(x,0)=sin (pi x).

Please see the accompanying PDF for more details about the code and its implementation.

***EDIT (17-Nov-2019):*** 
There was an error in the way the matrix-vector multiplication was written out in Eq.(13) [even though the code was correct], so this has been corrected in the equations and a new file has been re-uploaded, so as to avoid confusion. 
Thanks are due to Dr.Sainadh Chamarthi for not just pointing out the mistakes, but also for giving detailed explanation :)
