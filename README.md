Nonlinear Heat Transfer Simulation:

Description:
Simulates nonlinear heat transfer in a 1D rod using Implicit Euler time discretization and Finite Difference Method (FDM) in space. Newton’s method solves the nonlinear system arising from temperature-dependent thermal conductivity. The project visualizes both the temporal evolution and final temperature distribution for multiple nonlinearity parameters (χ).

Features:
Handles temperature-dependent thermal conductivity k(u) = ko e^χ&^u
Implements Implicit Euler method for stable time stepping
Uses Newton-Raphson iterations to solve nonlinear systems
Generates 3D mesh plots of temperature vs. space and time
Generates 2D plots of final temperature distribution along the rod

Requirements:
MATLAB R2018a or later (tested on R2023a)
No additional toolboxes required

Parameters:
rho, Cp : density and specific heat

kappa0 : baseline thermal conductivity

chi_values : array of nonlinearity parameters

a, b : spatial domain bounds

N : number of spatial grid points

tEnd, M : simulation end time and number of time steps

alpha, beta : Dirichlet boundary conditions at x=a and x=b

tol, maxIt : Newton-Raphson stopping tolerance and max iterations

Usage:
Open MATLAB.
Load the file temperature_plots_same_as_paper.m.

Run the function temperature_plots_same_as_paper():
View generated plots for temperature evolution and final distributions.

References:
Implicit Euler Time Discretization and Finite Difference Method with Newton’s Method in Nonlinear Heat Transfer Modeling, arXiv:1811.06337v1
