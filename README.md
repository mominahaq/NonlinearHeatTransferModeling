# Nonlinear Heat Transfer Simulation

## Description
This project simulates nonlinear heat transfer in a one-dimensional rod using the Implicit Euler method for time discretization combined with the Finite Difference Method (FDM) in space. The nonlinear system arising from temperature-dependent thermal conductivity is solved using Newton's method.

The simulation visualizes both the temporal evolution and the final temperature distribution along the rod for various nonlinearity parameters (χ).

## Features
- Models temperature-dependent thermal conductivity:  
  \( k(u) = k_0 e^{\chi u} \)
- Uses the Implicit Euler method for stable and accurate time stepping
- Solves nonlinear equations via Newton-Raphson iterations
- Produces 3D mesh plots of temperature as a function of space and time
- Generates 2D plots showing the final temperature distribution along the rod

## Requirements
- MATLAB R2018a or later (tested on R2023a)
- No additional toolboxes required

## Parameters
- `rho`, `Cp` : Density and specific heat capacity
- `kappa0` : Baseline thermal conductivity coefficient
- `chi_values` : Array of nonlinearity parameters χ
- `a`, `b` : Spatial domain boundaries (start and end points of the rod)
- `N` : Number of spatial grid points
- `tEnd`, `M` : Simulation end time and number of time steps
- `alpha`, `beta` : Dirichlet boundary conditions at \( x = a \) and \( x = b \)
- `tol`, `maxIt` : Newton-Raphson solver tolerance and maximum iteration count

## Usage
1. Open MATLAB.
2. Load the file `temperature_plots_same_as_paper.m`.
3. Run the function:
   ```matlab
   temperature_plots_same_as_paper()
