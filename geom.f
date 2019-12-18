	    MODULE geom
! (* 	----------------------------------------- *)
! (* 	This module contains the geometric information *)
! (* 	for the simulation. *)
! 
! (* 	Jorge Escobar-Vargas *)
! (* 	Cornell University *)
! (* 	February 2008 *)
! (* 	----------------------------------------- *)
	    USE scrotum
	    implicit none
	    save
	
! (* 	1. Array with coordinates of grid points *)
	    real, dimension(ngp,2) :: cgp
	
! (* 	2. Array with subdomain corner (nodal) points *)
	    integer, dimension(numsub,4) :: scp
	
! (* 	3. Array with type of Boundary condition  *)
	    integer, dimension(4) :: cond
	
! (* 	4. Array with the value of BC at the physical boundaries *)
	    real, dimension(4) :: val
	
! (* 	5. Integer with variable that sets periodic informatios *) (APR)
!	Suponiendo que: 0 es no periodico, 1 periodico en x, 2 periodico
! 	en z y 3 periodico en x y z
	    integer :: periodic = 1

! (*    6. Arreglo donde seran almacenadas las velocidades tomadas del
!	modelo numerico (APR - Febrero 2017)
	    real, allocatable, dimension(:,:) :: velocidades
	
	    END MODULE geom
