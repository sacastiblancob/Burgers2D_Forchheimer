	subroutine patching(B,B1,delta)
! 	----------------------------------------------------
! 	In this subroutine penalized patching conditions are 
! 	imposed at the subdomain interfaces. No storage
! 	
! 	VARIABLES
! 	taul -> Tau coefficient of left subdomain
! 	taur -> Tau coefficient of right subdomain
! 	taut -> Tau coefficient of top subdomain
! 	taub -> Tau coefficient of bottom subdomain
! 	lx -> Length of the left subdomain of the interface
! 	lx1 -> Length of the right subdomain of the interface
! 	lz -> Length of the bottom subdomain of the interface
! 	lz1 -> Length of the top subdomain of the interface
! 	
! 	Jorge Escobar-Vargas
! 	Cornell University - CEE
! 	February 2008
! 	----------------------------------------------------
!       MODULES
	USE scrotum
	USE geom
	USE legendre
	
	implicit none

! 	Dummy Variables	
	real, intent(in) :: delta
	real, dimension(nsg), intent(in) :: B
	real, dimension(nsg), intent(inout) :: B1
	
! 	Local Variables
	integer :: i,j,k,l
	integer :: r1,r11,r2,r22
	real :: lz,lz1,lx,lx1
	real :: alpha,beta,tau,omega,kappa,fac
	real :: taul,taur,taut,taub,caray,fbd
	integer :: temp

!	Se agrega la variable temp como esta en el codigo pasado por Escobar

! 	Penalty parameters
	alpha = 1.0
	beta = 1.0
	fac = 5.0e+2 !0.5 ! This factor is very important
	omega = 2. / (pd * (pd + 1.))
	kappa = omega * alpha / beta
	caray = (1. / (omega * delta * beta))
	tau = caray * (delta + (2. * kappa) -
     >	(2. * sqrt(kappa ** 2. + (delta * kappa))))
     
! 	1a)Continuity in the function u for vertical patching conditions
	if (nsubx /= 1) then
	do k = 0,nsubz-1
	  do j = 1,nsubx-1
	    lx = abs(cgp(scp(j,2),1) - cgp(scp(j,1),1))
	    lx1= abs(cgp(scp(j+1,2),1) - cgp(scp(j+1,1),1))
	    taul = tau * (2./lx) ! OJO I do not understand why this is not square
	    taur = tau * (2./lx1)! OJO I do not understand why this is not square
	    do i = 1,n
	      r1 = ns*((k*nsubx)+j-1)+(i*n)
	      r2 = ns*((k*nsubx)+j)+((i-1)*n)+1
! 	      Diagonal components of right condition in left subdomain
	      B1(r1) = B1(r1) - (fac * taul * alpha * B(r1))
! 	      Diagonal components of left condition in right subdomain
	      B1(r2) = B1(r2) - (fac * taur * alpha * B(r2))
!      	      Contribution of right subdomain to the condition of the left one
	      B1(r1) = B1(r1) + (fac * taul * alpha * B(r2))
! 	      Contribution of left subdomain to the condition of the right one
	      B1(r2) = B1(r2) + (fac * taur * alpha * B(r1))
	    enddo
	  enddo
	enddo
	endif
	
! 	1b)Continuity in the function u for horizontal patching conditions
	if (nsubz /= 1) then
	do k = 1,nsubz-1
!         El valor absoluto calculado en la siguiente linea se bajo APR
	  lz = cgp(scp((k-1)*nsubx+1,2),2) - cgp(scp((k-1)*nsubx+1,3),2)
          lz = abs(lz)
	  lz1 = abs(cgp(scp(k*nsubx+1,2),2) - cgp(scp(k*nsubx+1,3),2))
	  do j = 1,nsubx
	    taub = tau * (2./lz) ! OJO I do not understand why this is not square
	    taut = tau * (2./lz1)! OJO I do not understand why this is not square
	    do i = 1,n
	      r1 = ns*(((k-1)*nsubx)+j)-n+i
	      r2 = ns*((k*nsubx)+j-1)+i
! 	      Diagonal components of top condition in bottom subdomain
	      B1(r1) = B1(r1) - (fac * taub * alpha * B(r1))
! 	      Diagonal components of bottom condition in top subdomain
	      B1(r2) = B1(r2) - (fac * taut * alpha * B(r2))
! 	      Contribution of top subdomain to the condition of the bottom one
	      B1(r1) = B1(r1) + (fac * taub * alpha * B(r2))
! 	      Contribution of bottom subdomain to the condition of the top one
	      B1(r2) = B1(r2) + (fac * taut * alpha * B(r1))
	    enddo
	  enddo
	enddo
	endif
	
! 	2a) Continuity in the first derivative for vertical patching conditions
! 	In this case just the derivative with respect to x is taken into account
	if (nsubx /= 1) then
	fbd  = fac * beta * delta
	do k = 0,nsubz-1
	  do j = 1,nsubx-1
	    lx = abs(cgp(scp(j,2),1) - cgp(scp(j,1),1))
	    lx1 = abs(cgp(scp(j+1,2),1) - cgp(scp(j+1,1),1))
	    taul = tau * (2./lx)
	    taur = tau * (2./lx1)
	    do i = 1,n
	      do l = 0,n-1
	        r1  = ns*((k*nsubx)+j-1)+(i*n)
		r11 = ns*((k*nsubx)+j-1)+(i*n)-n+l+1
		r2  = ns*((k*nsubx)+j)+((i-1)*n)+1
		r22 = ns*((k*nsubx)+j)+((i-1)*n)+1+l
! 		Derivative of the right boundary of left subdomain
	        B1(r1) = B1(r1) - (fbd*taul*d(n,l+1)*(2./lx)*B(r11))
! 		Contribution of right subdomain to left one
		B1(r1) = B1(r1) + (fbd*taul*d(1,l+1)*(2./lx1)*B(r22))
! 		Contribution of left subdomain to right one
		B1(r2) = B1(r2) - (fbd*taur*d(n,l+1)*(2./lx)*B(r11))
! 		Derivative of the left boundary of right subdomain
		B1(r2) = B1(r2) + (fbd*taur*d(1,l+1)*(2./lx1)*B(r22))
	      enddo
	    enddo
	  enddo
	enddo
	endif
	
! 	2b) Continuity in the first derivative for horizontal patching conditions
	if (nsubz /= 1) then
	fbd  = fac * beta * delta
	do k = 1,nsubz-1
!         Valor absoluto movido para que no haya conflictos APR
	  lz = cgp(scp((k-1)*nsubx+1,2),2) - cgp(scp((k-1)*nsubx+1,3),2)
          lz = abs(lz)
	  lz1 = abs(cgp(scp(k*nsubx+1,2),2) - cgp(scp(k*nsubx+1,3),2))
	  do j = 1,nsubx
	    taub = tau * (2./lz)
	    taut = tau * (2./lz1)
	    do i = 1,n
	      do l = 0,n-1
	        r1  = ns*(((k-1)*nsubx)+j)-n+i
		r11 = ns*(((k-1)*nsubx)+j-1)+(l*n)+i
		r2  = ns*((k*nsubx)+j-1)+i
		r22 = ns*((k*nsubx)+j-1)+(l*n)+i
! 		Derivative of the top boundary of bottom subdomain
	        B1(r1) = B1(r1) - (fbd*taub*d(n,l+1)*(2./lz)*B(r11))
! 		Contribution of upper subd. to lower one
		B1(r1) = B1(r1) + (fbd*taub*d(1,l+1)*(2./lz1)*B(r22))
! 		Contribution of bottom subd. to upper one
		B1(r2) = B1(r2) - (fbd*taut*d(n,l+1)*(2./lz)*B(r11))
! 		Derivative of the bottom boundary of top subdomain
		B1(r2) = B1(r2) + (fbd*taut*d(1,l+1)*(2./lz1)*B(r22))
	      enddo
	    enddo
	  enddo
	enddo
	endif
	
! 	Desde aca se pega lo correspondiente a resolver NS con condiciones
! 	de contorno periodicas (APR)

!	************************************************************************
!	COMMENTS ON THIS PART CAN BE VERY CONFUSING
!	Implementing periodic boundary conditions for vertical boundaries
	if (periodic == 1 .or. periodic == 3) then
!	 Continuity in the function
	 do k = 0,nsubz-1
	  lx  = abs(cgp(scp((k+1)*nsubx,2),1) - cgp(scp((k+1)*nsubx,1),1))! Right Boundary
	  lx1 = abs(cgp(scp((k*nsubx)+1,2),1) - cgp(scp((k*nsubx)+1,1),1))! Left Boundary
	  taul = tau * (2./lx) 
	  taur = tau * (2./lx1)
	  do i = 1,n
	   r1 = ns*(((k+1)*nsubx)-1)+(i*n) ! Right Boundary
	   r2 = ns*(k*nsubx)+((i-1)*n)+1   ! Left Boundary

! 	   Diagonal components of right condition in left subdomain
	   B1(r1) = B1(r1) - (fac * taul * alpha * B(r1))
! 	   Diagonal components of left condition in right subdomain
	   B1(r2) = B1(r2) - (fac * taur * alpha * B(r2))
!      	   Contribution of right subdomain to the condition of the left one
	   B1(r1) = B1(r1) + (fac * taul * alpha * B(r2))
! 	   Contribution of left subdomain to the condition of the right one
	   B1(r2) = B1(r2) + (fac * taur * alpha * B(r1))
	  enddo
	 enddo
	 
! 	2a) Continuity in the first derivative for vertical patching conditions
! 	In this case just the derivative with respect to x is taken into account
	 fbd  = fac * beta * delta
	 do k = 0,nsubz-1
	  lx  = abs(cgp(scp((k+1)*nsubx,2),1) - cgp(scp((k+1)*nsubx,1),1))! Right Boundary
	  lx1 = abs(cgp(scp((k*nsubx)+1,2),1) - cgp(scp((k*nsubx)+1,1),1))! Left Boundary
	  taul = tau * (2./lx)
	  taur = tau * (2./lx1)
	  do i = 1,n
	   do l = 0,n-1
	    r1  = ns*(((k+1)*nsubx)-1)+(i*n) ! Right Boundary
	    r11 = ns*(((k+1)*nsubx)-1)+(i*n) - n + l + 1
	    r2  = ns*(k*nsubx)+((i-1)*n)+1   ! Left Boundary
	    r22 = ns*(k*nsubx)+((i-1)*n)+1 + l
	    
! 	    Derivative of the right boundary of left subdomain
	    B1(r1) = B1(r1) - (fbd*taul*d(n,l+1)*(2./lx)*B(r11))
! 	    Contribution of right subdomain to left one
	    B1(r1) = B1(r1) + (fbd*taul*d(1,l+1)*(2./lx1)*B(r22))
! 	    Contribution of left subdomain to right one
	    B1(r2) = B1(r2) - (fbd*taur*d(n,l+1)*(2./lx)*B(r11))
! 	    Derivative of the left boundary of right subdomain
	    B1(r2) = B1(r2) + (fbd*taur*d(1,l+1)*(2./lx1)*B(r22))
	   enddo
	  enddo
	 enddo
	
	endif	
!	************************************************************************

!	Implementing periodic BC for horizontal boundaries
	if (periodic == 2 .or. periodic == 3) then
!	 Continuity in the function
	  do k = 0,nsubx-1
	     temp = numsub-nsubx+1+k
	     lz  = abs(cgp(scp(temp,2),2) - cgp(scp(temp,3),2)) ! Top BC
	     lz1 = abs(cgp(scp(k+1,2),2)  - cgp(scp(k+1,3),2))  ! Bottom BC
	     taub = tau * (2./lz)
	     taut = tau * (2./lz1)
	     do i = 1,n
	       r1 = nsg-((nsubx-1)*ns)+(ns*k)-n+i !Top Boundary
	       r2 = (ns*k)+i ! Bottom Boundary
	       
! 	       Diagonal components of top condition in bottom subdomain
	       B1(r1) = B1(r1) - (fac * taub * alpha * B(r1))
! 	       Diagonal components of bottom condition in top subdomain
	       B1(r2) = B1(r2) - (fac * taut * alpha * B(r2))
! 	       Contribution of top subdomain to the condition of the bottom one
	       B1(r1) = B1(r1) + (fac * taub * alpha * B(r2))
! 	       Contribution of bottom subdomain to the condition of the top one
	       B1(r2) = B1(r2) + (fac * taut * alpha * B(r1))

	     enddo
	  enddo

!	 Continuity in the derivative
	  fbd  = fac * beta * delta
	  do k = 0,nsubx-1
	     temp = numsub-nsubx+1+k
	     lz  = abs(cgp(scp(temp,2),2) - cgp(scp(temp,3),2)) ! Top BC
	     lz1 = abs(cgp(scp(k+1,2),2)  - cgp(scp(k+1,3),2))  ! Bottom BC
	     taub = tau * (2./lz)
	     taut = tau * (2./lz1)
	     do i = 1,n
	       do l = 0,n-1

		r1  = nsg-((nsubx-1)*ns)+(ns*k)-n+i !Top Boundary
!		r11 = nsg-((nsubx-1)*ns)+(ns*k)-n+i-ns+((l+1)*n)
		r11 = nsg-((nsubx-k)*ns)+(n*l)+i
		r2  = (ns*k)+i ! Bottom Boundary
		r22 = (ns*k)+i+(l*n)

! 		Derivative of the top boundary of bottom subdomain
	        B1(r1) = B1(r1) - (fbd*taub*d(n,l+1)*(2./lz)*B(r11))
! 		Contribution of upper subd. to lower one
		B1(r1) = B1(r1) + (fbd*taub*d(1,l+1)*(2./lz1)*B(r22))
! 		Contribution of bottom subd. to upper one
		B1(r2) = B1(r2) - (fbd*taut*d(n,l+1)*(2./lz)*B(r11))
! 		Derivative of the bottom boundary of top subdomain
		B1(r2) = B1(r2) + (fbd*taut*d(1,l+1)*(2./lz1)*B(r22))

	      enddo

	     enddo
	  enddo
	endif

	end subroutine patching
