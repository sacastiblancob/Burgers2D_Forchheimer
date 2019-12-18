	subroutine CFL(u,w,coux,couz)
! 	------------------------------------------------------
! 	This subroutine calculates the Courant-Friedrichs-Lewy
! 	condition (CFL) for 2D Non-linear advection problems
! 	
! 	Jorge Escobar-Vargas
! 	Cornell University - CEE
! 	March 2008
! 	------------------------------------------------------
! 	MODULES
	USE scrotum
	USE map
	USE aetas
	
	implicit none
	real,dimension(nsg),intent(in) :: u,w
	real, intent(out) :: coux,couz
	
	integer :: i,j,k,loc,r1,r2
	real :: cout,dx
	
	coux = -1.0
	couz = -1.0
	
! 	Calculating CFL number for interior points
	
	do k = 0,numsub-1
	 do i = 1,n-2
	  do j = 2,n-1
	   loc = j + (i*n) + (k*ns)

! 	   Horizontal direction
	   
	   if (u(loc) >= 0) then
	     dx = abs(cx(loc+1) - cx(loc))
	   else
	     dx = abs(cx(loc) - cx(loc-1))
	   endif
! 	   write(*,*)k,i,j, dx
	   cout = dT * u(loc) / dx
	   coux = max(cout,coux)

! 	   Vertical direction

	   if (w(loc) >= 0) then
	     dx = abs(cz(loc+n) - cz(loc))
	   else
	     dx = abs(cz(loc) - cz(loc-n))
	   endif
	   cout = dT * w(loc) / dx
	   couz = max(cout,couz)
	   
	  enddo
	 enddo
	enddo
	
! 	Calculating CFL for Physical Boundaries
! 	Bottom Boundary
	do k = 0,nsubx-1
	 do i = 1,n
	  loc = (ns*k)+i
	  if (w(loc) > 0.0) then
	    dx = abs(cz(loc+n) - cz(loc))
	  endif
	  cout = dT * w(loc) / dx
	  couz = max(cout,couz)
	 enddo
	enddo  
	  
! 	Top Boundary  
	do k = 0,nsubx-1
	 do i=1,n
	   loc = nsg-(ns*k)+1-i
	   if (w(loc) < 0.0) then
	    dx = abs(cz(loc-n) - cz(loc))
	   endif
	   cout = dT * w(loc) / dx
	   couz = max(cout,couz)
	 enddo
	enddo
	
! 	Right Boudary
	do k = 1,nsubz
	  do i = 1,n
	    loc = ns*((k*nsubx)-1)+(i*n)
	    if (u(loc) < 0.0) then
	      dx = abs(cx(loc) - cx(loc-1))
	    endif
	    cout = dT * u(loc) / dx
	    coux = max(cout,coux)
	  enddo
	enddo
	
! 	Left Boundary     
	do k = 0,nsubz-1
	  do i = 0,n-1
	    loc = ns*(k*nsubx)+(i*n)+1
	    if (u(loc) > 0.0) then
	      dx = abs(cx(loc+1) - cx(loc))
	    endif
	    cout = dT * u(loc) / dx
	    coux = max(cout,coux)
	  enddo
	enddo
	
! 	Calculating CFL for interfaces
! 	Vertical Interfaces
	if (nsubx /= 1) then
	do k = 0,nsubz-1
	  do j = 1,nsubx-1
	    do i = 1,n
	      r1 = ns*((k*nsubx)+j-1)+(i*n)
	      r2 = ns*((k*nsubx)+j)+((i-1)*n)+1
	      if (u(r1) < 0.0) then
	        dx = abs(cx(r1) - cx(r1-1)) 
	      else
	        dx = abs(cx(r2+1) - cx(r2))
	      endif
	      cout = dT * u(r1) / dx
	      coux = max(cout,coux)
	      cout = dT * u(r2) / dx
	      coux = max(cout,coux)
	    enddo
	  enddo
	enddo
	endif
	
! 	Horizontal Interfaces
	if (nsubz /= 1) then
	do k = 1,nsubz-1
	  do j = 1,nsubx
	    do i = 1,n
	      r1 = ns*(((k-1)*nsubx)+j)-n+i
	      r2 = ns*((k*nsubx)+j-1)+i
	      if (w(r1) < 0.0) then
	        dx = abs(cz(r1) - cz(r1-n)) 
	      else
	        dx = abs(cz(r2+n) - cz(r2))
	      endif
	      cout = dT * w(r1) / dx
	      couz = max(cout,couz)
	      cout = dT * w(r2) / dx
	      couz = max(cout,couz)
	    enddo
	  enddo
	enddo
	endif
	
	end subroutine CFL