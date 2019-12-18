	subroutine error(t,u,w,ue,ve,Linfu,Linfw,du,dw)
! 	****************************************
! 	This subroutine calculates the Linf norm
! 	for the 2D Burgers equation
! 	
! 	Jorge Escobar-Vargas
! 	Cornell University - CEE
! 	June 2008
! 	****************************************
! 	MODULES
	
	USE scrotum
	USE aetas
	
	implicit none
	
	real, dimension(nsg), intent(in) :: u,w,ue,ve
	integer, intent(in) :: t
	real, intent(out) :: Linfu, Linfw
	
	real,dimension(nsg), intent(out) :: du,dw
	integer :: i
	
	do i=1,nsg
	 du(i) = abs(ue(i)-u(i))
	 dw(i) = abs(ve(i)-w(i))
	enddo
	
	Linfu = maxval(du)
	Linfw = maxval(dw)
	
! 	write(*,*)'After ',t*dT,' units of time'
! 	write(*,*)'Linf in u velocity is = ', Linfu
! 	write(*,*)'Linf in w velocity is = ', Linfw
	
	end subroutine error