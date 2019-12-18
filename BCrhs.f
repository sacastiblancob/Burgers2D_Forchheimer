	subroutine BCrhs(delta,BG,bv,idc)
! 	-----------------------------------------------------------
! 	In this subroutine penalized Dirichlet and Neumann boundary
! 	conditions are imposed for Helmholtz (Viscous) equation.
! 	In momentum equations
! 
! 	Jorge Escobar-Vargas
! 	Cornell University
! 	February 2008
! 	-----------------------------------------------------------
! 	MODULES
	USE scrotum
	USE geom
	USE legendre
	USE mound

	implicit none
	
! 	Dummy Variables
	integer, intent(in) :: idc
	real, intent(in) :: delta
	real,dimension(nsg),intent(inout) :: BG
	real,dimension(nsg),intent(in) :: bv
	
! 	Local Variables
	integer :: k, i, temp, temp1, nr, SS
	real :: lz,lx
	real :: tau, omega, alpha, beta, temp2, temp3

!	Verificando valores en el fondo del dominio computacional APR
	if (idc == 1) then
		temp3 = val(1)
	else if (idc == 2) then
		temp3 = val(1)
	endif
	
! 	Square domain (Global)

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
 	       BG(temp) = BG(temp) - (tau * temp3)
!	       BG(temp) = BG(temp) - (tau * bv(temp))
	     enddo
	   enddo
	elseif (cond(1) == 3) then ! Neumann with respect to z on the bottom
	   do k=0,nsubx-1
	     lz=abs(cgp(scp(k+1,2),2)-cgp(scp(k+1,3),2))
	     beta=1.
	     omega=2./(pd*(pd+1.))
	     tau=1.*(2./lz)/(beta*omega)
	     do i=1,n
	       nr = (ns*k)+i
	       BG(nr) = BG(nr) - (tau * delta * val(1))
	     enddo
	   enddo
	else
	  write(*,*) 'Something wrong in the bottom boundary conditions'
	endif

! 	Conditions for the top
	if (cond(3) == 1) then ! Dirichlet

		SS = size(velocidades(:,1))

	    do k = 0, nsubx-1
	     lz = abs(cgp(scp(numsub-k,2),2) - cgp(scp(numsub-k,3),2))
	     alpha = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = (delta*(2./lz)**2.) / (alpha*omega**2.)
	     do i=1,n
		   temp = nsg - (ns * k) + 1 - i
		   temp1 = SS
!	       temp1 = SS - (k + 1) * n + i
! 	       BG(temp) = BG(temp) - (tau * val(3))
!	       BG(temp) = BG(temp) - (tau * velocidades(temp1,idc))
	       temp2 = tau * velocidades(temp1, idc)
		   BG(temp) = BG(temp) - temp2
		   SS = SS - 1
	     enddo
	   enddo
	elseif (cond(3) == 3) then ! Neumann with respect to z (PUEDE MEJORARSE)
	   do k = 0,nsubx-1
	     lz = abs(cgp(scp(numsub-k,2),2) - cgp(scp(numsub-k,3),2))
	     beta = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = 1.*(2./lz) / (beta*omega)
	     do i = 1,n
	       nr = nsg-(ns*k)-n+i
	       BG(nr) = BG(nr) - (tau * delta * val(3))
	     enddo
	   enddo
	else
	  write(*,*) 'Something wrong in the top boundary conditions'
	endif
	endif

!	Esto es lo que se puso, hay que ver de donde viene PERIODIC
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
! 	       BG(temp) = BG(temp) - (tau * val(2))
	       BG(temp) = BG(temp) - (tau * bv(temp))
	     enddo
	   enddo
	elseif (cond(2) == 2) then ! Neumann with respect to x
	   do k = 1,nsubz
	     lx = abs(cgp(scp(k*nsubx,2),1)-cgp(scp(k*nsubx,1),1))
	     beta = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = 1.*(2./lx) / (beta*omega)
	     do i = 1,n
	       nr = ns*((k*nsubx)-1)+(i*n)
	       BG(nr) = BG(nr) - (tau * delta * val(2))
	     enddo
	   enddo
	else
	  write(*,*) 'Something wrong in the right boundary conditions'
	endif

! 	Conditions for the left
	if (cond(4) == 1 .or. cond(4) == 3) then ! Dirichlet
	   do k = 0,nsubz-1
!            El valor absoluto se saco de esta linea y se calculo luego APR	     
             lx = cgp(scp((k*nsubx)+1,2),1) - cgp(scp((k*nsubx)+1,1),1)
             lx = abs(lx)
	     alpha = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = (delta*(2./lx)**2.) / (alpha*omega**2.)
	     do i = 0,n-1
	       temp = ns*(k*nsubx)+(i*n)+1
! 	       BG(temp) = BG(temp) - (tau * val(4))
	       BG(temp) = BG(temp) - (tau * bv(temp))
	     enddo
	   enddo
	elseif (cond(4) == 2) then ! Neumann with respect to x
	   do k = 0,nsubz-1
!            El valor absoluto se saco de esta lnea y se puso abajo	     
             lx = cgp(scp((k*nsubx)+1,2),1) - cgp(scp((k*nsubx)+1,1),1)
             lx = abs(lx)
	     beta = 1.
	     omega = 2. / (pd*(pd+1.))
	     tau = 1.*(2./lx) / (beta*omega)
	     do i = 0,n-1
	       nr = ns*(k*nsubx)+(i*n)+1
	       BG(nr) = BG(nr) - (tau * delta * val(4))
	     enddo
	   enddo
	else
	  write(*,*) 'Something wrong in the left boundary conditions'
	endif
	endif
	
	end subroutine BCrhs
