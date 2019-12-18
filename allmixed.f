	subroutine allmixed(t,delta,B,B1)
! 	-------------------------------------------------
! 	This subroutine perform the Matrix-Vector product 
! 	of the global matrix (no BC, no patch) within the 
! 	solver (GMRES) avoiding matrix storage.

! 	Jorge Escobar-Vargas 
! 	Cornell University
! 	February 2008
! 	-------------------------------------------------
! 	MODULES
	USE scrotum
	USE geom
	USE legendre

	implicit none
	
! 	Dummy Variables	
	integer, intent(in) :: t
	real, intent(in) :: delta
	real, dimension(nsg), intent(in) :: B
	real, dimension(nsg), intent(out) :: B1
	
! 	Local Variables
	integer :: i,j,k,p
	integer :: row,colx,colz
	real :: lx, lz, suma
	
	do p = 0,numsub-1
! 	 Calculating the width and height of each subdomain
	 lx = abs(cgp(scp(p+1,2),1) - cgp(scp(p+1,1),1))
	 lz = abs(cgp(scp(p+1,2),2) - cgp(scp(p+1,3),2))

! 	 Matrix-Vector multiplication for the main matrix
! 	 Spatial derivatives (x and z)
	 
	 do k = 0,n-1
	  do i = 0,n-1
	    row = (p*ns)+(k*n)+i+1
	    suma=0.0
	    do j = 0,n-1
	      colx = (p*ns) + (k*n) + j + 1
	      colz = (p*ns) + (j*n) + i + 1
! 	      Second derivative with respect to x
	      suma = suma + (d2(i+1,j+1)*(2./lx)**2.0 * B(colx) * delta)
! 	      Second derivative with respect to z
	      suma = suma + (d2(k+1,j+1)*(2./lz)**2.0 * B(colz) * delta)
	    enddo
	    B1(row) = suma
	  enddo
	 enddo
	enddo 

	do i = 1,nsg
	  B1(i) = B1(i) - B(i)
	enddo
	
	end subroutine allmixed