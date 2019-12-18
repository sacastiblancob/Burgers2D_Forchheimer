	subroutine ex2dbur(ue,ve)
! 	---------------------------------------------
! 	This subroutine calculates the exact solution
! 	for the Steady 2D Burgers' equation
! 	
! 	See Fletcher(1991) Computational Techniques
! 	for Fluid Dynamics. Springer. Page 361
! 	
! 	Jorge Escobar-Vargas
! 	Cornel University - CEE
! 	February 2008
! 	---------------------------------------------
! 	MODULES
	USE scrotum
	USE mound
	USE map
	
	implicit none
! 	Dummy variables
	real, dimension(nsg), intent(out) :: ue,ve

! 	CON ESTO IMPONGO UNA CONDICION INICIAL DE VELOCIDAD CONSTANTE EN EL 
!	DOMINIO APR (170228) 	
	integer :: i

	do i = 1,nsg
	  ue(i) = vi_x
	  ve(i) = vi_z
	enddo

! 	CON ESTO PRUEBO LA SOLUCION DE LA GAUSSIANA MOVIENDOSE POR EL DOMINIO EN 
!	EL SENTIDO HORIZONTAL APR (170228)
! 	Local Variables
!	real :: A,b0,c0,S,d0,x0,z0
!	integer :: i

!	A = 0.5
!	x0 = 20.
!	z0 = -10.
!	S = 5.

!!	write(*,*) 'ddd', A, x0, z0, S
!!	call sleep(4)

!	do i=1,nsg
!	 
!	 b0 = ((cx(i) - x0)**2.) / (2. * S ** 2.)
!	 c0 = ((cz(i) - z0)**2.) / (2. * S ** 2.)
!	 d0 = -(b0 + c0)
!	 ue(i) = A * exp(d0)	 
!	 ve(i) = 0

!!	 write(*,*) 'COND INIC.', ue(i), ve(i), b0, c0
!	 
!	enddo

	
	
! 	c5 = exp(lambda*(-1.0-Xo)) + exp(-lambda*(-1.0-Xo))
! 	c6 = a1 + a2*(-1.0)+a3*(3.141592654 / 150.0)+
!      >	     a4*(-1.0)*(3.141592654 / 150.0)
! 	c7 = cos(lambda * 3.141592654 / 150.0)
! 	c8 = sin(lambda * 3.141592654 / 150.0)
! 	c9 = a3 + a4*(-1.0)
! 	lambda = c1 * (c9 - (c2*c5*c8)) / (c6 + (a5*c5*c7))
! 	
! 	write(*,*) 'Answer is: ', lambda
! 	
! 	write(*,*) nu,maxval(ue),maxval(ve)
! 	stop
	
	end subroutine ex2dbur
