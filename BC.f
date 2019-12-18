	subroutine BC(B,B1,delta)
! 	----------------------------------------------
! 	In this subroutine Dirichlet and Neumann  
! 	conditions are imposed for Helmholtz equation.
! 	No Storage
! 
! 	Jorge Escobar-Vargas
! 	Cornell University - CEE
! 	February 2008
! 	----------------------------------------------
! 	MODULES
	USE scrotum
	USE geom
	USE legendre
	
	implicit none
	
! 	Dummy Variables
	real, intent(in) :: delta
	real, dimension(nsg), intent(in) :: B
	real, dimension(nsg), intent(inout) :: B1
	
! 	Local Variables
	integer :: k,i,j,temp,nr,nc
	real :: lz,lx,tbd
	real :: tau,omega,alpha,beta

	if (periodic == 0 .or. periodic == 1) then	
! 	Conditions for the bottom
	if (cond(1) == 1) then ! Dirichlet
	   do k = 0,nsubx-1
	     lz = abs(cgp(scp(k+1,2),2) - cgp(scp(k+1,3),2))
	     alpha = 1.
	     omega = 2./(pd*(pd+1.))
	     tau = (delta*(2./lz)**2.) / (alpha*omega**2.)
	     do i = 1,n
	       temp = (ns*k)+i
	       B1(temp) = B1(temp) - (tau * alpha * B(temp))
	     enddo
	   enddo
	elseif (cond(1) == 3) then ! Neumann with respect to z on the bottom
	   do k=0,nsubx-1
	     lz=abs(cgp(scp(k+1,2),2)-cgp(scp(k+1,3),2))
	     beta=1.
	     omega=2./(pd*(pd+1.))
	     tau=1.*(2./lz)/(beta*omega)
	     tbd = tau * beta * delta
	     do i=1,n
	       do j=0,n-1
	         nr = (ns*k)+i
		 nc = (k*ns)+(j*n)+i
		 B1(nr) = B1(nr) + (tbd * d(1,j+1) * (2./lz) * B(nc))
	       enddo
	     enddo
	   enddo
	else
	  write(*,*) 'Something wrong in the bottom boundary conditions'
	endif

! 	Conditions for the top
	if (cond(3) == 1) then ! Dirichlet
! 	Aca meter un loop que sea de nsubx a 1
	   do k = 0,nsubx-1
	     lz = abs(cgp(scp(numsub-k,2),2) - cgp(scp(numsub-k,3),2))
	     alpha = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = (delta*(2./lz)**2.) / (alpha*omega**2.)
	     do i=1,n
	       temp = nsg-(ns*k)+1-i
	       B1(temp) = B1(temp) - (tau * alpha * B(temp))
	     enddo
	   enddo
	elseif (cond(3) == 3) then ! Neumann with respect to z (PUEDE MEJORARSE)
	   do k = 0,nsubx-1
	     lz = abs(cgp(scp(numsub-k,2),2) - cgp(scp(numsub-k,3),2))
	     beta = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = 1.*(2./lz) / (beta*omega)
	     tbd = tau * beta * delta
	     do i = 1,n
	       do j = n,1,-1
	         nr = nsg-(ns*k)+1-i
		 nc = nsg-(k*ns)+1-i-ns+(j*n)
		 B1(nr) = B1(nr) - (tbd * d(n,j) * (2./lz) * B(nc))
	       enddo
	     enddo
	   enddo
	else
	  write(*,*) 'Something wrong in the top boundary conditions'
	endif
	endif
	
!	Este if viene con la varibale PERIODIC
	if (periodic == 0 .or. periodic == 2) then
! 	Conditions for the right
	if (cond(2) == 1 .or. cond(2) == 3) then ! Dirichlet
	   do k = 1,nsubz
	     lx = abs(cgp(scp(k*nsubx,2),1) - cgp(scp(k*nsubx,1),1))
	     alpha = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = (delta*(2./lx)**2.) / (alpha*omega**2.)
	     do i = 1,n
	       temp = ns*((k*nsubx)-1)+(i*n)
	       B1(temp) = B1(temp) - (tau * alpha * B(temp))
	     enddo
	   enddo
	elseif (cond(2) == 2) then ! Neumann with respect to x
	   do k = 1,nsubz
	     lx = abs(cgp(scp(k*nsubx,2),1)-cgp(scp(k*nsubx,1),1))
	     beta = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = 1.*(2./lx) / (beta*omega)
	     tbd = tau * beta * delta
	     do i = 1,n
	       do j = 1,n
	         nr = ns*((k*nsubx)-1)+(i*n)
		 nc = ns*((k*nsubx)-1)+(i*n)-n+j
		 B1(nr) = B1(nr) - (tbd * d(n,j) * (2./lx) * B(nc))
	       enddo
	     enddo
	   enddo
	else
	  write(*,*) 'Something wrong in the right boundary conditions'
	endif

! 	Conditions for the left
	if (cond(4) == 1 .or. cond(4) == 3) then ! Dirichlet
	   do k = 0,nsubz-1
!            Quitar el valor ansoluto y calcularlo en la sguiente linea APR
	     lx = cgp(scp((k*nsubx)+1,2),1) - cgp(scp((k*nsubx)+1,1),1)
             lx = abs(lx)
	     alpha = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = (delta*(2./lx)**2.) / (alpha*omega**2.)
	     do i = 0,n-1
	       temp = ns*(k*nsubx)+(i*n)+1
	       B1(temp) = B1(temp) - (tau * alpha * B(temp))
	     enddo
	   enddo
	elseif (cond(4) == 2) then ! Neumann with respect to x
	   do k = 0,nsubz-1
!            quitar valor absoluto de esta linea y calcular en la siguiente APR
	     lx = cgp(scp((k*nsubx)+1,2),1) - cgp(scp((k*nsubx)+1,1),1)
             lx = abs(lx)
	     beta = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = 1.*(2./lx) / (beta*omega)
	     tbd = tau * beta * delta
	     do i = 0,n-1
	       do j = 1,n
	         nr = ns*(k*nsubx)+(i*n)+1
		 nc = ns*(k*nsubx)+(i*n)+j
		 B1(nr) = B1(nr) + (tbd * d(1,j) * (2./lx) * B(nc))
	       enddo
	     enddo
	   enddo
	else
	  write(*,*) 'Something wrong in the left boundary conditions'
	endif
	endif
	
	end subroutine BC
