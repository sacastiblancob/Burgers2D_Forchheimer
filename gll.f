! ****************************************************************************
! 	This file calculates the values of the Gauss-Lobatto-Legendre points
! 	within the interval [-1,1]
! 	This code is almost exactly the same than the code presented in Appendix C
! 	of the book "Spectral Methods in Fluid Dynamics" Canuto et.al 1982

!	JUST TO REMEMBER:   'n' is the number of points in each direction!!!


c**********************************************************************
      subroutine jacobl(n,alpha,beta,xcol,ndim)
!     subroutine jacobl(n,alpha,beta,xjac)  ! Original subroutine has this statement
c**********************************************************************
c
c  Computes the gauss-lobatto collocation points for Jacobi polynomials
c
c   n:              Degree of approximation (order of polynomials)
c   alpha:          Parameter in Jacobi weight
c   beta:           Parameter in Jacobi weight
c   xjac:           Output array that contains the roots from largest to smallest
c   xcol:           Output array with the collocation points
c   ndim:           Dimension of the array "xcol"
c
c
c  for Chebyshev-Gauss-Lobatto points use alpha=-0.5 and beta=-0.5
c  for Legendre-Gauss-Lobatto points use           0             0
c
      dimension xjac(1000),xcol(0:ndim)
!       real xjac(1)
      common /jacpar/alp,bet,rv
      data kstop/10/
      data eps/1.0e-12/
c
      alp = alpha
      bet = beta
      rv = 1. + alp
      np = n + 1
      call jacobf(np,pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1,1.0)
      call jacobf(np,pnp1m,pdnp1m,pnm,pdnm,pnm1m,pdnm1,-1.0)
      det = pnp*pnm1m - pnm*pnm1p
      rp = -pnp1p
      rm = -pnp1m
      a = (rp*pnm1m - rm*pnm1p)/det
      b = (rm*pnp   - rp*pnm)/det
      xjac(1) = 1.0
      nh = (n+1)/2
      dth = 3.14159265/(2*n+1)
      cd = cos(2.*dth)
      sd = sin(2.*dth)
      cs = cos(dth)
      ss = sin(dth)
c
      do 39 j=2,nh
       x = cs
       do 29 k=1,kstop
        call jacobf(np,pnp1,pdnp1,pn,pdn,pnm1,pdnm1,x)
        poly = pnp1 + a*pn + b*pnm1
        pder = pdnp1 + a*pdn + b*pdnm1
        recsum = 0.0
        jm = j - 1
        do 27 i=1,jm
          recsum = recsum + 1.0/(x-xjac(i))
27      continue
28      continue
        delx = -poly/(pder-recsum*poly)
        x = x +delx
        if(abs(delx) .lt. eps) go to 30
29      continue
30      continue
        xjac(j) = x
        cssave = cs*cd - ss*sd
        ss = cs*sd + ss*cd
        cs = cssave
39      continue
        xjac(np) = -1.0
        npp = n + 2
        do 49 i=2,nh
          xjac(npp-i) = -xjac(i)
49      continue
! 	if (n.ne. 2*(n/2)) go to 56
! 	xjac(nh+1)=0
! 	This part of the code was added by ???
        if(n.ne. 2*(n/2)) then 
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
         go to 56
        else
        xjac(nh+1) = 0.0
        do k=0,n
         kk = n - k + 1
         xcol(k) = xjac(kk)
         enddo
         return
        endif
56      return
        end
c
c**********************************************************************
        subroutine jacobf(n,poly,pder,polym1,pderm1,polym2,pderm2,x)
c**********************************************************************
! 	Computes the Jacobi polynomial (in this case Legendre) and its derivative
! 	of degree n at x
! 	
! 	n -> 	Polynomial degree
! 	poly->	Magnitude of Jacobi (Legendre) polynomial evaluated in "x "
! 	pder->	Magnitude of the derivative of Jacobi polynomial in "x"
! 	polym1->Polynomial of degree n-1
! 	pderm1->Derivative of the polynomial of degree n-1
! 	polym2->Polynomial of degree n-2
! 	pderm2->Derivative of the polynomial of degree n-2
! 	x ->	Point in the space
! 	
! 	This subroutine is an unabridged version of the original one presented in
! 	Canuto et. al. (1982)
	
	common /jacpar/ alp,bet,rv
        apb = alp + bet
        poly = 1.0
        pder = 0.0
        if(n.eq.0) return
        polylst = poly
        pderlst = pder
        poly = rv*x
        pder = rv
        if(n.eq.1) return
        do 19 k=2,n
          a1 = 2.*k*(k+apb)*(2.*k+apb-2.)
          a2 = (2.*k+apb-1.)*(alp**2-bet**2)
          b3 = (2.*k+apb-2.)
          a3 = b3*(b3+1.)*(b3+2.)
          a4 = 2.*(K+alp-1.)*(k+bet-1.)*(2.*k+apb)
          polyn = ((a2+a3*x)*poly - a4*polylst) / a1
          pdern = ((a2+a3*x)*pder - a4*pderlst + a3*poly) / a1
          psave = polylst
          pdsave = pderlst
          polylst = poly
          poly = polyn
          pderlst = pder
          pder = pdern
19      continue
        polym1 = polylst
        pderm1 = pderlst
        polym2 = psave
        pderm2 = pdsave
        return
        end

