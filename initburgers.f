	subroutine initburgers(ue,ve,u,w)
! 	**********************************************
! 	This subroutine generates the initial velocity
! 	fields for the 2D Burgers equation
! 	
! 	See Fletcher (1991) Computational Techniques
! 	for Fluid Dynamics I, page 360
! 	
! 	Jorge Escobar-Vargas
! 	Cornell University - CEE
! 	June 2008
! 	**********************************************
! 	MODULES
	USE scrotum
	USE map
        USE mound
	
	implicit none
	
	real, dimension(nsg), intent(in) :: ue,ve
	real, dimension(nsg), intent(out) :: u,w
	
	integer :: temp,i,k,temp1,l,temp2,j
	real :: m
	
	u = vi_x
	w = vi_z
	
! 	Initial conditions for the right and left boundary
	do k = 1,nsubz
	 do j = 1,nsubx
	  do i = 1,n
	   
	   temp = ns*((k*nsubx)-1)+(i*n) ! Right
	   u(temp) = ue(temp)
	   w(temp) = ve(temp)
	   
	   temp1 = ns*((k-1)*nsubx)+((i-1)*n)+1 ! Left
	   u(temp1) = ue(temp1)
	   w(temp1) = ve(temp1)
	    
	   do l = 1,n
	    temp2 = ns*((k-1)*nsubx)+(ns*(j-1))+(n*(i-1))+l
	    
	    if (temp2 /= temp .and. temp2 /= temp1 ) then
	     u(temp2) = 2. + ((-0.04-2.)/2.)*(cx(temp2)+1.)
	     m = (0.0 - w(temp1))/2.0
	     w(temp2) = m*(cx(temp2)-1.)
	    endif
	    
	   enddo
	  enddo
	  
	  
	 enddo
	enddo
	
! 	Initial conditions for the bottom boundary
	
	do k = 0,nsubx-1
	 do i = 1,n
	  temp = (ns*k)+i
!	SE CAMBIAN CONDICIONES INICIALES EN EL BOTTOM APR (170228)
!	  u(temp) = 0.
!	  w(temp) = 0.
	  u(temp) = ue(temp)
	  w(temp) = ve(temp)
	 enddo
	enddo
	
! 	Initial conditions for the top boundary

	do k = 0,nsubx-1
	 do i=1,n
	  temp = nsg-(ns*k)+1-i
	  u(temp) = ue(temp)
	  w(temp) = ve(temp)
	 enddo
	enddo
	
	end subroutine initburgers
