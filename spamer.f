	subroutine spamer(u,w,dudx,dudz,dwdx,dwdz,stx,stz)
! 	-------------------------------------------------------------
! 	This subroutine calculates the spatial term, which is used in
! 	the explicit scheme (AB) for Advection time advancement

! 	Jorge Escobar-Vargas 
! 	Cornell University
! 	February 2008
! 	-------------------------------------------------------------
! 	MODULES
	USE scrotum
	USE legendre
	
	implicit none
! 	Dummy Variables
	real,dimension(nsg),intent(in) :: u,w
	real,dimension(nsg),intent(in) :: dudx,dudz
	real,dimension(nsg),intent(in) :: dwdx,dwdz
	real,dimension(nsg),intent(out) :: stx,stz
! 	Local variables
	integer :: i
	real :: advx, advz
	real,allocatable,dimension(:) :: skewsym
	real,allocatable,dimension(:) :: temp1
	real,allocatable,dimension(:) :: temp2
!	Variables for reading Forchheimer parameters
	real, dimension(5) :: parameters
	real :: alfau, alfaw, betau, betaw

	open(unit=99, file='parameters.txt', status='old', action='read')
	read(99,*) parameters
	close(99)

	alfau = parameters(2)
	alfaw = parameters(3)
	betau = parameters(4)
	betaw = parameters(5)
	
!	write(*,*) alfau
!	write(*,*) betau

	if (advec == 1) then
	
! 	Primitive Variables form
	 write(*,*) 'PRIMITIVE VARIABLES'
	 do i = 1,nsg
	  stx(i) = (u(i)*dudx(i)) + (w(i)*dudz(i)) + alfau*u(i) + betau*u(i)*abs(u(i))
	  stz(i) = (u(i)*dwdx(i)) + (w(i)*dwdz(i)) + alfaw*w(i) + betaw*w(i)*abs(w(i))
!       Las siguientes dos lineas son una prueba no mas (APR)
!	  stx(i) = 3.0*dudx(i) 
!	  stz(i) = 3.0*dwdx(i) 
	 enddo
	
	elseif (advec == 2) then
	
! 	 Skew-symmetric form of the non-linear term

	 allocate(skewsym(nsg),temp1(nsg),temp2(nsg))
	
	 do i=1,nsg
	  skewsym(i) = 0.5 * u(i) * u(i)
	 enddo
	 call diffx(skewsym,stx) ! Here stx = 0.5 * duudx
	
	 do i=1,nsg
	  skewsym(i) = 0.5 * u(i) * w(i)
	 enddo
	 call diffx(skewsym,temp1)! Here temp1 = 0.5 * duwdx
	 call diffz(skewsym,temp2)! Here temp2 = 0.5 * duwdz
	
	 do i=1,nsg
	  skewsym(i) = 0.5 * w(i) * w(i)
	 enddo
	 call diffz(skewsym,stz) ! Here stz = 0.5 * dwwdz
	
	
	 do i=1,nsg
	
	  advx = 0.5 * ((u(i)*dudx(i)) + (w(i)*dudz(i)))
	  stx(i) = stx(i) + temp2(i) + advx
	 
	  advz = 0.5 * ((u(i)*dwdx(i)) + (w(i)*dwdz(i)))
	  stz(i) = stz(i) + temp1(i) + advz
	 
	 enddo

	 deallocate(skewsym,temp1,temp2)

	elseif (advec == 3) then 

	  ! Divergence form
	  allocate(skewsym(nsg),temp1(nsg),temp2(nsg))
	
	 do i=1,nsg
	  skewsym(i) = u(i) * u(i)
	 enddo
	 call diffx(skewsym,stx) ! Here stx = duudx
	
	 do i=1,nsg
	  skewsym(i) = u(i) * w(i)
	 enddo

	 call diffx(skewsym,temp1)! Here temp1 = duwdx
	 call diffz(skewsym,temp2)! Here temp2 = duwdz
	
	 do i=1,nsg
	  skewsym(i) = w(i) * w(i)
	 enddo
	 call diffz(skewsym,stz) ! Here stz = dwwdz
	
	
	 do i=1,nsg
	
	  stx(i) = stx(i) + temp2(i)
	 
	  stz(i) = stz(i) + temp1(i)
	 
	 enddo

	 deallocate(skewsym,temp1,temp2)
	
	endif
	
	end subroutine spamer
