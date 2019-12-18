	subroutine interavg2d(t,rho)
! 	************************************
! 	Strong adaptive interfacial average
! 	See Diamessis (2005) 298-322 JCP
! 	
! 	Jorge Escobar-Vargas
! 	Cornell University - CEE
! 	June 2008
! 	************************************
! 	MODULES
	USE scrotum
	USE aetas
	
	
	implicit none
	
	integer, intent(in) :: t
	real, parameter :: Cave = 0.01
	real, dimension(nsg), intent(inout) :: rho
	
	real :: C0,c1,c2,c3
	integer :: r1, r2, i, j, k
	integer :: cor1,cor2,cor3,cor4
	real,dimension(6) :: C
	integer :: countV, countVa
	integer :: countH, countHa
	integer :: countC, countCa
	
	real :: bu1,bu2
	
! 	Averaging in Vertical Interfaces
	
	if (nsubx /= 1) then
	countV = 0
	countVa = 0
	do k = 0,nsubz-1
	  do j = 1,nsubx-1
	    do i = 1,n
	      
	      countV = countV + 1
	    
	      r1 = ns*((k*nsubx)+j-1)+(i*n) ! -> ||
	      r2 = ns*((k*nsubx)+j)+((i-1)*n)+1 ! || <-
	      
	      c1 = abs(rho(r2)-rho(r1))
	      c2 = abs(rho(r2)+rho(r1))
	      
	      C0 = c1 / c2
	      
	      if (C0 > Cave) then
	       
	       countVa = countVa + 1
	      
	       c3 = 0.5 * (rho(r1-1)+rho(r2+1))
	       
	       rho(r1) = c3
	       rho(r2) = c3
	       
	      endif
	    enddo
	  enddo
	enddo
	endif
	
! 	Averaging Horizontal interfaces
	
	if (nsubz /= 1) then
	countH = 0
	countHa = 0
	do k = 1,nsubz-1
	  do j = 1,nsubx
	    do i = 1,n
	      
	      countH = countH + 1
	    
	      r1 = ns*(((k-1)*nsubx)+j)-n+i ! ^ =
	      r2 = ns*((k*nsubx)+j-1)+i ! v =
	      
	      c1 = abs(rho(r2)-rho(r1))
	      c2 = abs(rho(r2)+rho(r1))
	      
	      C0 = c1 / c2
	      
	      if (C0 > Cave) then
	      
	       countHa = countHa + 1

	       c3 = 0.5 * (rho(r1-n)+rho(r2+n))
	       
	       rho(r1) = c3
	       rho(r2) = c3
	       
	      endif
	    enddo
	  enddo
	enddo
	endif
	
! 	Averaging for corners
	
	if (nsubx /= 1 .and. nsubz/=1) then
	 countC = 0
	 countCa = 0
	 do k = 1,nsubz-1
	  do j = 1,nsubx-1
	   
	   countC = countC + 1
	  
	   cor1 = (ns*(k-1)*nsubx) + (j*ns)
	   cor2 = (ns*(k-1)*nsubx) + ((j+1)*ns) -n + 1
	   cor3 = (k*nsubx*ns) + ((j-1)*ns) + n
	   cor4 = (k*nsubx*ns) + (j*ns) + 1
	   
	   C(1) = abs(rho(cor1)-rho(cor4)) / abs(rho(cor1)+rho(cor4))
	   C(2) = abs(rho(cor2)-rho(cor3)) / abs(rho(cor2)+rho(cor3))
	   C(3) = abs(rho(cor1)-rho(cor3)) / abs(rho(cor1)+rho(cor3))
	   C(4) = abs(rho(cor2)-rho(cor4)) / abs(rho(cor2)+rho(cor4))
	   C(5) = abs(rho(cor1)-rho(cor2)) / abs(rho(cor1)+rho(cor2))
	   C(6) = abs(rho(cor3)-rho(cor4)) / abs(rho(cor3)+rho(cor4))
	   
	   C0 = maxval(C)
	   
	   if (C0 > Cave) then
	   
	       countCa = countCa + 1

	       c3 =      (rho(cor1-1)+rho(cor1-n))
	       c3 = c3 + (rho(cor2+1)+rho(cor2-n))
	       c3 = c3 + (rho(cor3-1)+rho(cor3+n))
	       c3 = c3 + (rho(cor4+1)+rho(cor2+n))
	       c3 = c3 / 8.0
	       
	       rho(cor1) = c3
	       rho(cor2) = c3
	       rho(cor3) = c3
	       rho(cor4) = c3
	       
	   endif
	   
	  enddo
	 enddo
	endif
	
! 	Statistics for the interfacial and corners average
	bu1 = (countVa+countHa)*100/(countV+countH)
	bu2 = countCa*100/countC
	
	if (bu1>0. .or. bu2>0.) then
	 write(*,*)
	 write(*,*)'t = ', t
! 	 write(*,*)'Total numbero of interfacial points = ',countV+countH
!	 write(*,*)'Averaging applied to ', bu1 ,'% of the points'
! 	 write(*,*)'Total number of corners = ',countC
!	 write(*,*)'Averaging applied to ',bu2 ,'% of corners'
	endif
	end subroutine interavg2d
