# SVD_2020
The MATLAB code is a finite difference approximation for the Fisher-Stefan model in n-dimensions to explore the spreading-vanishing dichotomy (SVD) This MATLAB file takes the discretized equations from Simpson (2020) and integrates forward in time using a standard central difference scheme, the resulting nonlinear systems of algebraic equations are solved using Newton-Raphson iteration and the resulting linear systems are solving using the Thomas algorithm.

The MAPLE worksheet contains the exact leading eigenvalue solution to the linearised partial differential equation, verification that the solution solves the governing equation, visualisation of the solution, and solution of the mass balance equation to give $L_{c}$. 
