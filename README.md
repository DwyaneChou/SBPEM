# SBPEM
A Spherical Barotropic Primitive Equation Model

This model had been written in Matlab code.

Main Characteristics
Equations        : IAP transformation
L_operator       : Inversly Symmetric Operator
B_operator       : Harmoniously Diffusive Operator
Time Integration : Variable time step square conservation integral scheme

There is a Rossby-Haurwitz wave as initial condition for test case, also you can use a 500hPa u,v,Z to run this model.
u : u component wind speed in m/s
v : u component wind speed in m/s
Z : Geopotential Height in m^2/s^2

The output file name is "output.nc", you can plot the fields rapidly by using Panoply,ncview,ncbrowse and so on.