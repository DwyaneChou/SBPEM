# SBPEM
# A Spherical Barotropic Primitive Equation Model

This model had been written in Matlab code.

Main Characteristics

Mesh             : Arakawa C grid

Equations        : Shallow water equations with IAP transformation.

L_operator       : Antisymmetry Operator

B_operator       : Harmoniously Diffusive Operator

Time Integration : 2nd order conservation integral scheme with variable time step

Integral Scheme  :

1. HDO(Harmoniously Diffusive Operator)

2. A improved 4th order Runge-Kuuta

3. Predict-Correct

4. Leap-frog

Split Scheme     : CSP2(2nd order Conservative Splitting Pattern)


There is a Rossby-Haurwitz wave as initial condition for test-case, also you can use the 500hPa u,v,Z to run this model.

u : u component wind speed in m/s

v : v component wind speed in m/s

Z : Geopotential Height in m^2/s^2

h : sqrt(Z), the gravity wave phase speed

U : h\*u

V : h\*v

# File Intrduction

SBPEM.m is the main script for this model;

L_operator.m is the antisymmetry operator;

IAP.m is the IAP transformation program;

fast_pass.m and slow_pass.m are the fast-wave part and slow-wave part for split scheme;

genMesh.m is the mesh generator;

B_operator.m is the (HDO)Harmoniously Diffusive Operator;

HDO.m is the HDO integrate scheme;

LF.m is the Leap-Frog integrate scheme;

PC is the Predict-Correct integrate scheme;

RK4 is a improved Runge-Kutta scheme;

Haurwitz.m is the script to generate the Rossby-Haurwitz wave as the initial condition for testing the model;

split_integrator.m is the split integrate procedure;

spatial_discrete.m is the script which be used to choose wave-pass type;

inner_product.m is used to calculate the inner product;

output_netCDF.m is the I/O program.

The output file name is "output.nc", you can plot the fields rapidly by using Panoply, ncview, ncbrowse and so on.


Zhou Lilong

National Meteorological Center of CMA

College of Earth and Planetary Sciences, UCAS

Special thanks to my teacher Dong Li @IAP:https://github.com/dongli

# Update log
V3.2 :

Fixed a bug in L operator, now the model can run much more better than before.

V3.1 :

(1) Optimized code structure, integrated the variables into 3 structures: MESH, STATE and TEND.

(2) Improved the running speed.

(3) Fixed several bugs.

V3.0 :

(1) Added integral scheme RK4, Leap-frog, Predict-Correct.

(2) Added 2nd order conservation split scheme.

(3) Modified the grid. Now, the mass grid is on the full grid, and v is on the half grid.
