	subroutine BDAB(t,rho,st,rhom1,rhom2,stm1,stm2)
! 	----------------------------------------------
! 	In this subroutine is implemented the
! 	Adams-Bashforth method to advance in time the
! 	Advective component of momentum equations

! 	See Karniadakis(1991) JCP 97,414-443
! 	See Peyret(2002) Spectral Methods for 
! 	    Incompressible Viscous Flow. Springer (page 131)
! 	See Ascher(1995) SIAM J.Num.Anal. Vol 32 No 3.797-823
	   
! 	Jorge Escobar-Vargas
! 	Cornell University
! 	February 2008
! 	----------------------------------------------
! 	MODULES
	USE scrotum
	USE aetas
	
	implicit none
	
! 	Dummy variables
	integer, intent(in) :: t
	real, dimension(nsg), intent(in) :: st
	real, dimension(nsg), intent(inout) :: rho
	real, dimension(nsg), intent(inout) :: rhom1, rhom2
	real, dimension(nsg), intent(inout) :: stm1, stm2
! 	Local Variables
	real :: bd, bd1, bd2
	real :: ab, ab1, ab2
	real,allocatable,dimension(:) :: rhof

	allocate(rhof(nsg))
	
	 if (t == 1) then
	 
! 	   First fractional time step BDF1/AB1

	   rhof  = rho - (dT * st)
	   rhom1 = rho
	   rho = rhof
	   stm1 = st
	 
	 elseif (t == 2) then
	 
! 	   Second fractional Time step BDF2/AB2
	   
	   bd = 2.0
	   bd1 = -1.0/2.0
	   ab =  2.0
	   ab1 = -1.0
	   rhof  = (bd*rho) + (bd1*rhom1) - 
     >		   (dT * ((ab*st) + (ab1*stm1)))
	   rhom2 = rhom1
	   rhom1 = rho
	   rho = rhof
	   stm2 = stm1
	   stm1 = st
	 
	 else
	  
! 	   General fractional time-step BDF3/AB3
	   
	   bd = 3.0
	   bd1 = -3.0 / 2.0
	   bd2 = 1.0 / 3.0
	   ab = 3.0
	   ab1 = -3.0
	   ab2 = 1.0
	   rhof  = (bd*rho) + (bd1*rhom1) + (bd2*rhom2) -
     >	           (dT * ((ab*st) + (ab1*stm1) + (ab2*stm2)))
	   rhom2 = rhom1
	   rhom1 = rho
	   rho = rhof
	   stm2 = stm1
	   stm1 = st
	 
	 endif
	
	 deallocate(rhof)
	 
	end subroutine BDAB