	subroutine filtering(n,numsub,ns,nsg,rho,F)
	
! 	This subroutine performs the spectral filtering
! 	to the numerical solution.
! 	See Blackburn, JCP 186 (2003) 610 - 629 for details
! 	Jorge Escobar - Cornell University
! 	June 2007
	
	implicit none
	
	integer, intent(in) :: n, numsub, ns, nsg
	real, dimension(n,n), intent(in) :: F
	real, dimension(nsg), intent(inout) :: rho
	
! 	Local Variables
	
	real, dimension(nsg) :: rhoN
	real, dimension(n,n) :: rhoT, rhoF1, rhoF2
	integer :: i,j,k,c
	
	rhoN = 0.
	do k = 0, numsub - 1
	 rhoT = 0.0
	 rhoF1 = 0.0

! 	 Extracting the subdomain information from global solution
	 do i = 0, n - 1
	  do j = 1, n
	   c = (ns * k) + (i * n) + j
	   rhoT(i + 1, j) = rho(c)
	  enddo
	 enddo
! 	 Now, filtering the solution at the subdomain 
! 	 rhoT=transpose(rhoT)
	 rhoF1 = matmul(F, rhoT)
	 rhoF2 = matmul(rhoF1, transpose(F))
! 	 rhoF1=transpose(rhoF1)
! 	 Returning the filtered solution to the global system
	 do i = 0,n - 1
	  do j = 1, n
	   c = (ns * k) + (i * n) + j
	   rhoN(c) = rhoN(c) + rhoF2(i + 1, j)
	  enddo
	 enddo
	
	enddo
	rho = rhoN
	
	
	end subroutine filtering
