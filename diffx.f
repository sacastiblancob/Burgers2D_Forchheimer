	subroutine diffx(rho,dx)
! 	-----------------------------------------------
! 	This subroutine calculates the X - derivatives
! 	for the skew-symmetric advective term of the 
! 	Shallow Water Equations

! 	Jorge Escobar-Vargas
! 	Cornell University
! 	February - 2008
! 	-----------------------------------------------
! 	MODULES
	USE scrotum
	USE geom
	USE legendre

	implicit none
! 	Dummy Variables
	real, dimension(nsg), intent(in) :: rho
	real, dimension(nsg), intent(out) :: dx
! 	Local Variables
	real, dimension(n) :: R
	integer :: p, i, j, k
	real :: lx
	
	do k = 0,numsub-1
	 lx = abs(cgp(scp(k+1,2),1) - cgp(scp(k+1,1),1))
	 do j = 0,n-1
	  do i = 1,n
	   p = (k*ns) + (j*n) + i
	   R(i) = rho(p)
	  enddo
	  R = matmul(d,R)
	  R = R * (2.0 / lx)
	  do i = 1,n
	   p = (k*ns) + (j*n) + i
	   dx(p) = R(i)
	  enddo
	 enddo
	enddo
	
	end subroutine diffx