c**********************************************************************
      subroutine quad(n,x,w,ndim)
c**********************************************************************
!       include 'dim.h'
!       include 'comlegendre.h'
      parameter (nn=1000)
      dimension alpol(0:n,0:n),dalpol(0:n,0:n)
      dimension x(0:ndim),w(0:ndim),alp1(0:nn),al1(0:nn)
c
c  **PD-modify**: Also store values of 0 to ndim Legendre
c  polynomials at all collocation points.
c  determine the Gauss Quadrature weighting factors
c
      small = 1.0e-30
      do k=0,n
       xc = x(k)
       call legen(al1,alp1,n,xc,nn)  
c-Calculate weighting factor
       w(k) = 2. / 
     &         ( n*(n+1)*al1(n)*al1(n) + small )
c-Store values for j-th Legendre polynomial at point xc=x(k)
c-Store values for derivative of j-th order L. polynomial at
c-collocation point, as well.
       do j=0,nzm
         alpol(j,k) = al1(j)
         dalpol(j,k) = alp1(j)
       enddo

      enddo

      return
      end

c**********************************************************************
      subroutine legen(al,alp,n,xc,ndim)
c**********************************************************************
C----------------------------------------------------------------------
C-Calculates values of all Legendre polynomials (and immediately highest
C-one in hierarchy) at a given collocation point xc.
C-Same operation for derivatives of Legendre polynomials.
      dimension al(0:ndim),alp(0:ndim)
c
      al(0) = 1.
      al(1) = xc
      alp(0) = 0
      alp(1) = 1.
c
      do k=1,n
        kp = k + 1
        km = k - 1
        al(kp) = (2.*k+1.)*xc*al(k)/kp - k*al(km)/kp
      enddo
c
      do k=1,n
        kp = k + 1
        km = k - 1
        alp(kp) = (2.*k+1.)*(al(k)+xc*alp(k))/kp - 
     &            k*alp(km)/kp
      enddo
c
      return
      end
c
      