	subroutine diffz(rho,dz)
! 	-----------------------------------------------
! 	This subroutine calculates the Z - derivatives
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
	real, dimension(nsg), intent(out) :: dz
! 	Local Variables
	real, dimension(n) :: R
	integer :: p, i, j, k
	real :: lz
	
	
	do k = 0,numsub-1
	 lz = abs(cgp(scp(k+1,2),2) - cgp(scp(k+1,3),2))
	 do j = 1,n
	  do i = 0,n-1
	   p = (k*ns) + j + (i*n)
	   R(i+1) = rho(p)
	  enddo
	  R = matmul(d,R)
	  R = R * (2.0 / lz)
	  do i = 0,n-1
	   p = (k*ns) + j + (i*n)
	   dz(p) = R(i+1)
	  enddo
	 enddo
	enddo
	
	end subroutine diffz