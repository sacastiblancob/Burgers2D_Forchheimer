	    MODULE legendre
! 	-----------------------------------------------------
! 	This module defines the Legendre collocation points
! 	as well as the differentiation matrices associated to 
! 	those points.
! 	
! 	Jorge Escobar-Vargas
! 	Cornell University
! 	February 2008
! 	-----------------------------------------------------
	    USE scrotum
	    implicit none
	
! 	1. Array with Gauss-Lobatto-Legendre collocation points
	    real, allocatable, dimension(:) :: points
	
! 	2. 1st, 2nd and 3rd order differentiation matrices
	    real, allocatable, dimension(:,:) :: d, d2, d3

	    END MODULE legendre