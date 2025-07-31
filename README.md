# TransPore Version 2.0 2023

•	A FEM-FDM solver for two-phase, multicomponent transport in porous media with non-Newtonian (shear-thinning) effects for polymer component

GENERAL USAGE NOTES
=====================
•	Requirements: MATLAB 2017a or higher, Windows/Linux/MacOS
•	E-mail: daripa@math.tamu.edu
•	Copyright 2010-2018 TransPore developers and contributors. All rights reserved.
•	References: 
•	Daripa, P. & Mishra, R. (2023) Modelling shear thinning polymer flooding using a dynamic viscosity model arXiv preprint arXiv:2301.04290
•	Daripa, P. & Dutta, S. (2017) Modeling and simulation of surfactant–polymer flooding using a new hybrid method. J. Comput. Phys., 335, 249–282.
https://doi.org/10.1016/j.jcp.2017.01.038
•	Daripa, P. & Dutta, S. (2019) On the convergence analysis of a hybrid method for multicomponent transport in porous media, Appl. Numer. Math., 146 (2019), pp. 199-220. https://doi.org/10.1016/j.apnum.2019.07.009


•	Funding for this research was provided by:
•	Qatar National Research Fund (08-777-1-141)
•	National Science Foundation (DMS-1522782)

Primary SOURCE FILE
===================
Master_surf_grid.m
Variables and Data Structure:
nsim - number of different flooding simulations
sizeofgrid - Nx x Ny grid sizes for each simulation   
c0iter, g0iter - nsim arrays of concentrations of components 1 & 2 respectively in the injected fluid for each simulation
f - source term for the elliptic problem. Nx x Ny matrix with non-zero intensities at injection and production wells
KK - Nx x Ny matrix with absolute permeability values for the domain
UU, CC, GG - Nx x Ny matrices for space-time values of wetting phase saturation, components 1 & 2 concentrations respectively 
miuw, miuo - wetting and non-wetting fluid base viscosities respectively
swr0, sor0 - wetting and non-wetting phase initial residual saturations respectively
sigma - Nx x Ny matrix for interfacial tension values over the domain
miua - Nx x Ny matrix for aqueous phase saturation values over the domain
lamba_a, lambda_o - Nx x Ny matrices for wetting phase and non-wetting phase mobility values over the domain 
u,v - Nx x Ny matrices for x-direction and y-direction total velocity values over the domain. Note: These are obtained by solving the global pressure equation and are different from phase velocities

Secondary SOURCE FILES
=====================
KKdef() - function implementing different types of homogeneous, heterogeneous, stochastic, piecewise constant absolute permeability profiles
s0c0() - function implementing initial configurations and injection profiles for each simulation.
compvis(), compres(), compmob() - functions to update phase viscosities, phase residual saturations and phase mobility values with the evolution of state variables
setGrid(), setRightHand(), setA(), setB(), getu(), get_vn() - functions implementing various parts of the elliptic solver for the global pressure and total velocities
nmmoc_surf_mod_neumann() - function implementing the MMOC-FD procedure for solving the component transport equations 
  

Instructions to run the code with Non-Newtonian modeling

Follow these steps to run the non-Newtonian TransPore code:
1.	Unpack the code to your local directory
2.	The code is designed in such a way that minimum access is required to other sub-routines and most of the changes can be made from the master file (master_surf_grid.m).
3.	Different configurations can be set up based on the flags in the first 150 lines of the master_surf_grid.m which is the main file in the code.
4.	Although, there are lot of parameters that can be changed, the main ones that govern the simulations and are important are listed as :

```
nsim =1;
sog=29;
sizeofgrid = [sog sog sog sog];
c0iter = [0.001 0 0 0.001 ]; %0.0002 (300 wppm) 0.0006 (900 wpppm) 0.001 (1500 wpppm)
g0iter = [0 0 0 0]; %[0.01 0.01 0.01];
Here, if you wish you can run 4 or more simulations in one go. The nsim entity ensures the number of simulations to be executed (currently set to 1 which means only 1 simulation will run during code execution). For each of one of these simulations, different values of IPC (c0) and ISC (g0) can be assigned from the array c0iter and g0iter. The sog parameter is quite important as it sets the number of mesh cells in the domain. The number of cells in x and y are the same and is equal to the sog value. 
src = 120000
```

src is the IR which determines how fast the injected fluid will move through the domain. The values are selected in such a way that the generated shear is within the range of the empirical constants for the polymer. The range is between 10,000 to 120,000. 
```
    %% Rectilinear propagation
    f(:,1)=src;                           % intensity of injection well = src
    f(:,para.box.m+1)=-src;    % intensity of production well = -src

%     %% Quarter-five spot
%   f(1,1)=src;                           % intensity of injection well = src
%   f(para.box.n+1,para.box.m+1)=-src;    % intensity of production well = -src
```
This part of the code will select either rectilinear or quarter five spot geometry. Comment out the one that is not being used.

KKdef(counter,sizeofgrid,x,y,permeabilityFlag);

This defines the permeability matrix. Notice the permeability flag can be chosen according to the following:
```
    %%%%%%%   Defining permeability matrix    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flag values:  1 = homogeneous with magnitude 1000
    %               2 = continuous heterogeneous function
    %               3 = impermeable block inclusion at the center
    %               4 = impermeable block inclusion off the center
    %               5 = Upper Ness from SPE10 model sections
    %               6 = Tabert from SPE10 model sections
s0=0.79;  % initial residual water saturation inside the reservoir = 1-s0
1-s0 defines the initial water saturation. This means that the oil field is not fully composed of oil initially but a mixture of oil and water in 0.79 (oil) and 0.21 (water).

viscosityFlag = 3; %Flag to enable different viscosity models
    % 3 = Dynamic Viscosity Model (Non-Newtonian model)
    % 2 = Sourav's Model
    % 1 = Original model from JCP paper
    
polymerType = 0; %Flag to initialize with a particular polymer
    % 0 = Xanthane
    % 1 = Schizophyllan
```
These flags are essential as they indicate which viscosity model and polymer is chosen. To enable the non-Newtonian model set the viscosityFlag to 3.

Note: phi_test=get_phi_test(para); 

This calls a function get_phi_test which in turn calls z_func_test. To change the shape of the initial conditions and to move from rectilinear to quarter five spot go to the z_func_test and change the function based on the comments.

5.	Once you are confident with your setup, you can run the simulation. It will take around 2-3hrs to finish a simulation if the sog is 29. Increasing the value of sog increases the computational time exponentially.
6.	Main post-processing variables are COC, MFW, CC, UU, GG. Depending on what kind of results you want to analyze you can plot those variables. 
