!----------------------------------------------------------------------!
!This subroutine builds the vector product of ALHS*PHI for the SEM
!on Quadrilateral Elements for the Shallow Water Equations on the Sphere
!for the Geopotential Semi-Implicit.
!Written by Francis X. Giraldo on 8/99
!           Naval Research Laboratory
!           Global Modeling Section
!           Monterey, CA 93943-5502
!----------------------------------------------------------------------!
	subroutine lhs_gmres(alhs,x,t,delta)

	USE scrotum
	USE geom
	USE legendre
	
	implicit none
	
	integer, intent(in) :: t
	real, intent(in) :: delta
	real, dimension(nsg) :: alhs,x
	
	alhs = 0.0 ! Very important The new basis vector is created
	
! 	Global system of equations
	call allmixed(t,delta,x,alhs)

! 	Boundary conditions
	call BC(x,alhs,delta)
      
! 	Patching conditions
	call patching(x,alhs,delta)
	
	
	
	
      end subroutine lhs_gmres



