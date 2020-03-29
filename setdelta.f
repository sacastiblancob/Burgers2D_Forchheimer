	subroutine setdelta(t,rho,BG,delta)
! 	-----------------------------------------
! 	In this subroutine the gamma0 parameter
! 	from time-splitting is introduced into
! 	the Helmholtz (Viscous) term of the 
! 	momentum equation
! 	
! 	Jorge Escobar-Vargas
! 	Cornell University
! 	February 2008 
! 	-----------------------------------------
! 	MODULES
	USE scrotum
	USE mound
	USE aetas
	
	implicit none
	
! 	Dummy Variables
	integer,intent(in) :: t
	real, intent(out) :: delta
	real,dimension(nsg),intent(in) :: rho
	real,dimension(nsg),intent(out) :: BG
! 	Local Variables
	
	if (t == 1) then
	
	     delta = nu * dT / 1.0
	     
	elseif (t == 2) then
	
	     delta = nu * dT / (3.0 / 2.0)
	     
	else
	
	     delta = nu * dT / (11.0 / 6.0)
	     
	endif
	
	BG = -rho * delta / (nu * dT)
	
	end subroutine setdelta