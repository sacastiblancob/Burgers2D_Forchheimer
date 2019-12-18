	subroutine patchpen(u,w,rho,st)
! 	-------------------------------------------------
! 	This subroutine imposes patching conitions to the
! 	Advective terms of the momentum equations

! 	Jorge Escobar-Vargas 
! 	Cornell University
! 	February 2008

! 	LIST OF LOCAL VARIABLES
! 	taul -> Tau coefficient of left subdomain
! 	taur -> Tau coefficient of right subdomain
! 	taut -> Tau coefficient of top subdomain
! 	taub -> Tau coefficient of bottom subdomain
! 	lx -> Length of the left subdomain of the interface
! 	lx1 -> Length of the right subdomain of the interface
! 	lz -> Length of the bottom subdomain of the interface
! 	lz1 -> Length of the top subdomain of the interface
! 	mfl -> Mapping factor for LEFT subdomain
! 	mfr -> Mapping factor for RIGHT subdomain
! 	mfb -> Mapping factor for BOTTOM subdomain
! 	mft -> Mapping factor for TOP subdomain
! 	-------------------------------------------------

! 	MODULES
	USE scrotum
	USE geom
	
	implicit none
	
! 	Dummy Variables
	real,dimension(nsg),intent(in) :: rho,u,w
	real,dimension(nsg),intent(out) :: st !! OJO INOUT
! 	Local variables
	integer :: i,j,k,r1,r2
	real :: lx, lx1, lz, lz1, mfl, mfr, mfb, mft
	real :: alpha, gamma, omega
	real :: taul,taur,taut,taub,fac1,fac2
	
! 	Some parameters
	fac1 = 1.0
	fac2 = 1.0 
	omega = 2.0 / (pd * (pd + 1.0))

!  	Patching for vertical interfaces
	if (nsubx /= 1) then
	do k = 0,nsubz-1
	  do j = 1,nsubx-1
	    lx = abs(cgp(scp(j,2),1) - cgp(scp(j,1),1))
	    mfl =  2.0 / lx
	    lx1= abs(cgp(scp(j+1,2),1) - cgp(scp(j+1,1),1))
	    mfr = 2.0 / lx1
	    do i = 1,n
	      r1=ns*((k*nsubx)+j-1)+(i*n)     ! Left Subd
	      r2=ns*((k*nsubx)+j)+((i-1)*n)+1 ! Right Subd
	      if (u(r1) >= 0.0) then ! Apply penalty to Right side of the interface
		alpha = u(r1)
		taur = fac1 * mfr / (2.0 * omega)
		st(r2) = st(r2) + (taur * ((alpha*rho(r2)) - (alpha*rho(r1))))
	      else ! Apply penalty to LEFT side of the interface
	        gamma = abs(u(r1))
		taul = fac2 * mfl / (2.0*omega)
		st(r1) = st(r1) + (taul * ((gamma*rho(r1)) - (gamma*rho(r2))))
	      endif
	    enddo
	  enddo
	enddo
	endif

! 	Patching for horizontal interfaces
	if (nsubz /= 1) then
	do k = 1,nsubz-1  ! OJO, BAJAR ESTOS LZ PARA QUE CALCULE POR SUBDOMINIO.
	                  ! usar algo que relac numsub=nsubx*nsubz y hacerle
	  lz = cgp(scp((k-1)*nsubx+1,2),2) - cgp(scp((k-1)*nsubx+1,3),2)
!	  quite el valor absoluto en la linea 74, hice la diferencia y lo calcule despues APR	  
      lz = abs(lz)
	  mfb = 2.0 / lz
	  lz1 = abs(cgp(scp(k*nsubx+1,2),2) - cgp(scp(k*nsubx+1,3),2))
	  mft = 2.0 / lz1
	  do j = 1,nsubx
	    do i = 1,n
	      r1 = ns*(((k-1)*nsubx)+j)-n+i ! Bottom Subdomain
	      r2 = ns*((k*nsubx)+j-1)+i     ! Top Subdomain
	      if (w(r1) >= 0.0) then ! Apply penalty to Right side of the interface
		alpha = w(r1)
		taut = fac1 * mft / (2.0 * omega)
		st(r2) = st(r2) + (taut * ((alpha*rho(r2)) - (alpha*rho(r1))))
	      else ! Apply penalty to LEFT side of the interface
	        gamma = abs(w(r1))
		taub = fac2 * mfb / (2.0*omega)
		st(r1) = st(r1) + (taub * ((gamma*rho(r1)) - (gamma*rho(r2))))
	      endif
	    enddo
	  enddo
	enddo
	endif
	
	end subroutine patchpen
