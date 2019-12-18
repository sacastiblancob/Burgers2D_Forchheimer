	subroutine velocity(u,w)
! 	------------------------------------------
! 	This subroutine calculates the velocity
! 	fields u(x,z) and w(x,z) for the linear
! 	advection-diffusion process

! 	Jorge Escobar - Cornell University
! 	February 2008
! 	------------------------------------------

	USE scrotum
	USE map

	implicit none
	
! 	Dummy Arguments
	real, dimension(nsg), intent(out) :: u,w
	
! 	Local Variables
	integer :: i
	
! 	CHANGE EVERYTHING BELOW THIS POINT.

! 	For the circular bump
! 	do i = 1,nsg
! 	 u(i) = (cz(i)-0.0) ! ojo menos la mitad del dominio
! 	 w(i) = (-cx(i)+0.0) ! ojo mas la mitad del dominio 
! 	enddo
! 	Just advection in one direction
	do i = 1,nsg
	 u(i) = -0.5
	 w(i) = 0.0
	enddo
	
	end subroutine velocity