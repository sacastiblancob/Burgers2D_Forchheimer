# 2DBurgers - Forchheimer

This repository contains a Fortran program that solves the Burger's equation
in a rectangular two dimensional domain using the spectral multi domain penalty
method developed by Escobar at Cornell University as part of  PhD project. 

The program has been modified in order to gather velocity signals on the top of
the domain, periodic boundary conditions in the right and left boundaries and a 
Dirichlet BC in the bottom of the domain. 

There is no current reference available for the project. 

Antonio Preziosi Ribero

170420 - Initial push of the program with functions

The code has been modified in order to take into account the 4 Forchheimer terms
into de Burgers2D equation, in fact you Don't need to recompile the codes.
If you want yo change the viscosity or the Forchheimer parameters, you just should
modify the parameters.txt file.

1. Viscosity
2. Linear Forchheimer in U direction
3. Nonlinear Forchheimer in U direction
4. Linear Forchheimer in W direction
5. Nonlinear Forchheimer in W direction

200131 - Sergio Castiblanco



