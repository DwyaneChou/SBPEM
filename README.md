# SBPEM
A Spherical Barotropic Primitive Equation Model

This model had been written in Matlab code.

Main Characteristics

Equations        : IAP transformation

L_operator       : Antisymmetry Operator

B_operator       : Harmoniously Diffusive Operator

Time Integration : Square conservation integral scheme with variable time step

Integral Scheme  :

1. HDO(harmonious diffusion operator)

2. A improved 4th order Runge-Kuuta

3. Predict Correct

4. leap frog

Split Scheme     : CSP2(2nd order Conservative Splitting Pattern)


There is a Rossby-Haurwitz wave as initial condition for test-case, also you can use the 500hPa u,v,Z to run this model.

u : u component wind speed in m/s

v : v component wind speed in m/s

Z : Geopotential Height in m^2/s^2


The output file name is "output.nc", you can plot the fields rapidly by using Panoply,ncview,ncbrowse and so on.


Zhou Lilong

National Meteorological Center of CMA

College of Earth and Planetary Sciences, UCAS

Special thanks to my teacher Dong Li @IAP:https://github.com/dongli
