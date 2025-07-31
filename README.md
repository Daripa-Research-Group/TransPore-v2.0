# TransPore v2.0

A FEM-FDM solver for two-phase, multicomponent transport in porous media. This version includes non-Newtonian (shear-thinning) effects for the polymer component.

---

## 📋 General Information
* **Requirements**: MATLAB R2017a or higher on Windows, Linux, or macOS.
* **Contact**: For questions or support, please email `daripa@math.tamu.edu`.
* **Copyright**: © 2010-2025 TransPore developers and contributors. All rights reserved.

---

## 💰 Funding
This research was made possible by funding from:
* Qatar National Research Fund (08-777-1-141)
* National Science Foundation (DMS-1522782)

---

## 📚 Citing this Work
If you use this code in your research, please cite the relevant publications:

1.  **Daripa, P. & Mishra, R. (2023).** Modelling shear thinning polymer flooding using a dynamic viscosity model. *Physics of Fluids, 35(4)*. [https://doi.org/10.1063/5.0145061](https://doi.org/10.1063/5.0145061)
2.  **Daripa, P. & Mishra, R. (2023).** Modelling shear thinning polymer flooding using a dynamic viscosity model. *arXiv preprint arXiv:2301.04290*.
3.  **Daripa, P. & Dutta, S. (2017).** Modeling and simulation of surfactant–polymer flooding using a new hybrid method. *Journal of Computational Physics, 335*, 249–282. [https://doi.org/10.1016/j.jcp.2017.01.038](https://doi.org/10.1016/j.jcp.2017.01.038)
4.  **Daripa, P. & Dutta, S. (2019).** On the convergence analysis of a hybrid method for multicomponent transport in porous media. *Applied Numerical Mathematics, 146*, 199-220. [https://doi.org/10.1016/j.apnum.2019.07.009](https://doi.org/10.1016/j.apnum.2019.07.009)

---

## 📁 File Structure

### Primary Source File
The main file for setting up and running simulations. Most parameters are configured here.
* `Master_surf_grid.m`

### Secondary Source Files
These files contain helper functions called by the master script.
* `KKdef()`: Implements various permeability profiles (homogeneous, heterogeneous, etc.).
* `s0c0()`: Sets initial configurations and injection profiles.
* `compvis()`, `compres()`, `compmob()`: Update phase viscosities, residual saturations, and mobilities.
* `setGrid()`, `setRightHand()`, `setA()`, `setB()`, `getu()`, `get_vn()`: Implement parts of the elliptic solver for global pressure and total velocity.
* `nmmoc_surf_mod_neumann()`: Implements the MMOC-FD procedure for solving component transport equations.

---

## 🔬 Key Variables & Data Structures
The primary variables are stored as `Nx x Ny` matrices.
* `nsim`: Number of flooding simulations to run.
* `sizeofgrid`: An array specifying the `Nx x Ny` grid sizes for each simulation.
* `c0iter`, `g0iter`: Arrays of injected concentrations for components 1 & 2.
* `f`: Source term matrix for the elliptic problem, with non-zero values at injection/production wells.
* `KK`: Absolute permeability matrix for the domain.
* `UU`, `CC`, `GG`: Space-time values for wetting phase saturation (UU), component 1 concentration (CC), and component 2 concentration (GG).
* `miuw`, `miuo`: Base viscosities for the wetting and non-wetting fluids.
* `swr0`, `sor0`: Initial residual saturations for the wetting and non-wetting phases.
* `sigma`: Interfacial tension over the domain.
* `miua`: Aqueous phase viscosity over the domain.
* `lamba_a`, `lambda_o`: Wetting phase and non-wetting phase mobility values.
* `u`, `v`: Total velocity components in the x and y directions.

---

## 🚀 Running a Simulation

Follow these steps to configure and run the non-Newtonian TransPore code.

**1. Unpack the Code**: Extract the files to a local directory.

**2. Configure the Master File**: Open `master_surf_grid.m`. All primary simulation parameters can be adjusted within the first ~150 lines.

### Key Configuration Parameters

#### **Simulation Count & Grid Size**
Set the number of simulations (`nsim`) and the grid size (`sog`). The domain is a square grid of `sog` x `sog` cells.

```matlab
nsim = 1;
sog = 29; % Sets the number of mesh cells (e.g., 29x29)
sizeofgrid = [sog sog sog sog];
```

Injected Fluid Concentrations
Define the initial concentrations for each component for each simulation.
```matlab
% c0iter: Polymer concentration (e.g., 0.001 = 1500 wppm)
c0iter = [0.001 0 0 0.001];
% g0iter: Surfactant concentration
g0iter = [0 0 0 0];
```
Injection Rate & Well Geometry
src determines the injection rate (IR). The geometry (well placement) can be set to rectilinear or quarter-five spot.
```maltb
% Recommended range: 10,000 to 120,000
src = 120000;

% Choose a geometry by uncommenting the desired block
%% Rectilinear propagation
f(:,1) = src;          % Injection well
f(:,para.box.m+1) = -src; % Production well

%% Quarter-five spot
% f(1,1) = src;                                % Injection well
% f(para.box.n+1,para.box.m+1) = -src;         % Production well
```

Permeability Profile
Use the permeabilityFlag to select the domain's permeability profile.

```matlab
% 1 = Homogeneous (1000 mD)
% 2 = Continuous heterogeneous function
% 3 = Impermeable block at the center
% 4 = Impermeable block off-center
% 5 = Upper Ness (from SPE10 model)
% 6 = Tarbert (from SPE10 model)
permeabilityFlag = 3;
KKdef(counter, sizeofgrid, x, y, permeabilityFlag);
```

Initial Water Saturation
The value 1-s0 defines the initial water saturation in the reservoir.

```matlab
s0 = 0.79; % Initial oil saturation
% Corresponds to an initial water saturation of 1 - 0.79 = 0.21
```

Viscosity Model & Polymer Type
These flags are essential for enabling the non-Newtonian model.

```matlab
% Set viscosityFlag to 3 for the non-Newtonian model
viscosityFlag = 3;
% 3 = Dynamic Viscosity Model (Non-Newtonian)
% 2 = Sourav's Model
% 1 = Original JCP Model

% Set the polymer type
polymerType = 0;
% 0 = Xanthane
% 1 = Schizophyllan
```
Note: The initial saturation profile shape is defined in get_phi_test, which calls z_func_test. Modify z_func_test to change the initial conditions from rectilinear to quarter-five spot.

**3. Run Simulation**: Once your setup is configured, run master_surf_grid.m from MATLAB.

A simulation with sog=29 may take 2-3 hours.

Computational time increases exponentially with sog.

**4. Post-Processing**: Analyze the results by plotting the main output variables: COC (cost of chemicals), MFW (mass fraction of water), CC (component 1 conc.), UU (saturation), and GG (component 2 conc.).
