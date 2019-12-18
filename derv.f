! **********************************************************************
! **********************************************************************

! 	This code contains the routine employeed to calculate the first,
! 	second, and third derivative matrix at Gauss-Lobatto-Legendre
! 	grid points calculated in the routines jacobl and jacobf.
! 	
! 	The theory behind this code can be found in
! 	Costa & Don "On the computation of high order derivatives". 
! 	Applied numerical mathematics. 33 (2000) 151 - 159


!       call derv(nzm,zpts,d,d2,d3,nzm)
c**********************************************************************
      subroutine derv(nterm,x,d,d2,d3,ndim)
c*********************************************************************
! 	List of Variables:
! 	nterm ->	Polynomial degree
! 	x ->		Array with Gauss-Lobatto-Legendre Points
! 	d ->		First derivative matrix
! 	d2 ->		Second derivetive matrix
! 	d3 ->		Third derivative matrix
! 	ndim ->		Number of grid points in each direction
	
      parameter (nn=1000)
      dimension x(0:ndim),c(0:nn),
     >          d(0:ndim,0:ndim),d2(0:ndim,0:ndim),d3(0:ndim,0:ndim)
!      >          al1(0:nn),alp1(0:nn),
!      >          al2(0:nn),alp2(0:nn),c(0:nn)

C-Introduce methodology by Costa and Don.
!	see Eqn (4)
C-1) Must first compute coefficients C_i
      do 10 k=0,nterm
        prod = 1.
        xk = x(k)

        do 20 l=0,nterm
          xl = x(l)

          if (l.ne.k) then
            prod = prod*(xk-xl)
          else
            prod = prod
          endif

 20     continue
        c(k) = prod

 10   continue

C-This is the lengthy route, but no big deal.
C-FIRST DERIVATIVE   
C-Calculate off diagonal elements of D^1
!	see Eqn (6)
      do 30 k=0,nterm
        xk = x(k) 

        do 40 j=0,nterm
          xj = x(j)
          if (k.ne.j) then
            d(k,j) = c(k)/(c(j)*(xk-xj))
          else
            d(k,j) = 0.
          endif

 40     continue
 30   continue

C-Calculate diagonal elements now
C-Diagonal element is negative row sum of off diagonal
C-elements (See Costa & Don) Eqn(9)
      do 50 k=0,nterm
        sum = 0.
        do 60 j=0,nterm

          if (k.ne.j) then
            sum = sum + d(k,j)
          else
            sum = sum
          endif

 60     continue
        d(k,k) = -sum
c       if ((k.ne.0).and.(k.ne.nterm)) d(k,k) = 0.
 50   continue

C-SECOND  DERIVATIVE
C-Off-diagonal elements Eqn(13)
      m = 2
      do 70 k=0,nterm
        xk = x(k)

        do 80 j=0,nterm
          xj = x(j)
          if (k.ne.j) then
            d2(k,j) = float(m)*(d(k,k)*d(k,j) - d(k,j)/(xk-xj))
          else
            d2(k,j) = 0.
          endif

 80     continue
 70   continue      

C-Diagonal elements Eqn(9) again
       do 90 k=0,nterm
        sum = 0.
        do 100 j=0,nterm

          if (k.ne.j) then
            sum = sum + d2(k,j)
          else
            sum = sum
          endif

 100    continue
        d2(k,k) = -sum
 90   continue

C-THIRD  DERIVATIVE
C-Off-diagonal elements
      m = 3
      do 110 k=0,nterm
        xk = x(k)

        do 120 j=0,nterm
          xj = x(j)
          if (k.ne.j) then
            d3(k,j) = float(m)*(d2(k,k)*d(k,j) - d2(k,j)/(xk-xj))
          else
            d3(k,j) = 0.
          endif

 120    continue
 110  continue

C-Diagonal elements
       do 130 k=0,nterm
        sum = 0.
        do 140 j=0,nterm

          if (k.ne.j) then
            sum = sum + d3(k,j)
          else
            sum = sum
          endif

 140    continue
        d3(k,k) = -sum
 130  continue
      return
      end subroutine derv


