	MODULE scrotum
! (* 	------------------------------------------------------------ *)
! (* 	This MODULE contains the basic parameters for the simulation *)
! (* 	of the Shallow Water Equations via pseudospectral multidomain*)
! (* 	penalty method *)
! 
! (* 	Jorge Escobar-Vargas *)
! (* 	Cornell University *)
! (* 	February 2008 *)
! (* 	------------------------------------------------------------ *)
	implicit none
	save
	
! (* 	1. Number of collocation points in each direction per subdomain *)
	integer, parameter :: n = 5
	
! (* 	2. Number of subd in X direction  *)
	integer, parameter :: nsubx = 30
	
! (* 	3. Number of subd in Z direction  *)
	integer, parameter :: nsubz = 10
	
! (* 	4. Number of grid points *)
	integer, parameter :: ngp = (nsubx+1)*(nsubz+1)
	
! (* 	5. Total number of subdomains *)
	integer, parameter :: numsub = nsubx * nsubz
	
! (* 	6. Number of collocation points per subdomain (ns=n*n) *)
	integer, parameter :: ns = n * n
	
! (* 	7. Global number of collocation points (nsg=ns*numsub)  *)
	integer, parameter :: nsg = ns * numsub
	
! (* 	8. Polynomial degree *)
	integer, parameter :: pd = n - 1

! (*    9. Forma de tratar el termino advectivo de la ecuacion skew o
!	forma "natural" *)
	integer :: advec = 1
	
	END MODULE scrotum
