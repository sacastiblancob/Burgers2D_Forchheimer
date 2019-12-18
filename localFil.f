	subroutine localFil(n,ns,points,F,wg,pf)
! 	This subroutine generates the matrix for spectral
! 	filtering
! 	Jorge Escobar - Cornell University
! 	June 2007
	
	implicit none
	
	integer, intent(in) :: n,ns
	real, intent(in) :: pf
	real, dimension(n), intent(in) :: points, wg
	real, dimension(n,n), intent(out) :: F
	
! 	Local variables
	
	real, dimension(n,n) :: B,L,M,C,W
	real, dimension(n) :: P
	integer :: kp,k,km,i,j
	real :: alpha,pi,c1,c2
	real :: kf, nf
	real :: specfilter_bvfourier2d
	
	pi = acos(-1.0)
	nf = n

! 	B matrix -> Legendre polyn. evaluated in GLL points
	do j = 1, n
	 P(1) = 1.
	 P(2) = points(j)
	 do k = 2, n - 1
	  kp = k + 1
	  km = k - 1
	  P(kp)=(((2.0*(k-1))+1.)*P(2)*P(k)/(k+1-1))-((k-1)*P(km)/(k+1-1))
	 enddo
	 do i = 1, n
	  B(j, i) = P(i)
	 enddo
	enddo
	C = 0.0
	W = 0.0
	L = 0.0
	do i = 1, n - 1
	 C(i, i) = (i - 1.) + 0.5
	 W(i, i) = wg(i)
	enddo
	C(n,n) = n / 2.
	W(n,n) = wg(n)
	M = matmul(C, transpose(B))
	M = matmul(M, W)

! 	Exponential filter
	alpha = -log(1.0e-16) ! Machine precision according to Diamesis 2008
	write(*,*)'Alpha is: ', alpha
	write(*,*)'pf is: ', pf
	do k = 1,n
	  kf = k
	  L(k, k) = exp(-alpha * ((k - 1.0) / (n - 1.0)) ** pf)

! 	  L(k,k)=specfilter_bvfourier2d(kf,0.0,nf,pf)
	enddo

! 	c2 = sqrt(2.0)*(n-1)
! 	do i=1,n
! 	 do k=1,n
! 	  c1 = sqrt((i-1)**2.0 + (k-1)**2.0)
! 	  L(i,k) = exp(-alpha*(c1**2.0/c2)**pf)
! 	 enddo
! 	enddo
	
	F = matmul(B, L)
	F = matmul(F, M)
	end subroutine localFil
	
! Comente toda esta funcion por recomendacion de escobar (APR)

!***********************************************************************
!       function specfilter_bvfourier2d(ak,akcf,akmax,p)
!***********************************************************************
!-PD:5/29/02.
!-Implementation of Boyd-Vandeven spectral filter proposed by Boyd (1997)
!-and used by Levin et al. (JCP 137, 130-154 (1997)) in spectral element ocean circulation model.
!-Arguments:
!-a) k = Legendre wavenumber under consideration
!-b) kc = limit wavenumber above which filtering is allowable.
!-This is the lag. Choose kc=2/3*N as indicated in Taylor et al. (JCP 1997)
!-c) N = number of points in Legendre expansion
!-Filter order is set to p=2 to damp out noise at small scales.
! 
!      if ((ak.gt.akcf).and.(ak.lt.akmax)) then
! 
!       theta = (ak-akcf)/(akmax-akcf)
!!-Proper calculation of chi parameter.
!       omega = abs(theta) - 0.5
!       if (theta.eq.0.5) then
!         chi = 1.
!       else
!         c1 = 4.*omega**2.
!           c2 = -log(1.-c1)/c1
!           chi = sqrt(c2)
!       endif
!       arg = 2.*sqrt(p)*chi*omega
!       specfilter_bvfourier2d = 0.5*erfc(arg)
! 
!      else
!          if (ak.eq.akmax) then
!            specfilter_bvfourier2d = 0.
!          else
!          specfilter_bvfourier2d = 1.
!           endif
!      endif
! 
!      end function specfilter_bvfourier2d
