# 2DBurgers

This repository contains a Fortran program that solves the Burger's equation
in a rectangular two dimensional domain using the spectral multi domain penalty
method developed by Escobar at Cornell University as part of  PhD project. 

The program has been modified in order to gather velocity signals on the top of
the domain, periodic boundary conditions in the right and left boundaries and a 
Dirichlet BC in the bottom of the domain. 

There is no current reference available for the project. 

Antonio Preziosi Ribero

170420 - Initial push of the program with functions

With Forchheimer parameters taking into account, you DON'T need to compile everything everytime,
you should just change the parameters.txt file. The order of the parameters is:
Viscosity
Alfa u
Alfa w
Beta u
Beta w

Alfa is the first order Forchheimer, and Beta is the second order Frochheimer corresponding to
velocities U and W.


