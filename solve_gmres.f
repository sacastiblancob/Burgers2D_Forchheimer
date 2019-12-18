!-----------------------------------------------------------------------
!      GMRES Solver: 
!             Original Version by P.F. Fischer, Argonne National Lab
!             Modified Version by F.X. Giraldo, Naval Research Lab
! 	      See Saad - Iterative Methods for Sparse Linear Systems
! 	      Modified version by J.A. Escobar-Vargas for a specific case
! 	      NO matrix storage. See subroutine lhs_gmres
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
	subroutine solve_gmres(x,b,ts,delta,niter)      
      
	USE scrotum
	USE geom
	USE legendre
	
	implicit none
	
! 	Dummy Variables
	integer, intent(in) :: ts
	integer, intent(out) :: niter
	real, intent(in) :: delta
	real, dimension(nsg), intent(out) :: x
	real, dimension(nsg), intent(in) :: b
	
! 	Local Variables
	integer, parameter :: mmax = 600 ! Max number of iterations - Cambiar si necesario
	real, allocatable, dimension(:) :: adiag, w
	real, allocatable, dimension(:,:) :: v ! Basis Vectors
	real, allocatable, dimension(:,:) :: h ! Hessenberg Matrix
	real, allocatable, dimension(:) :: s,c,gamma
	
	real :: twonorm, glsc2
	logical ifdone
	
	real :: one,eps,solver_tol,gamma1,etol,etol2
	real :: t,den
	
	integer :: i,j,k,maxiter,kiter,kiter0,kiter1
	
	!set machine tolerances - aca se puede cambiar la tolerancia para esto
	solver_tol = 1.e-6
	one = 1.
	eps = 1.e-20
	if (one+eps == one) eps = 1.e-14
	if (one+eps == one) eps = 1.e-7
	if (one+eps == one) eps = 1.e-16
	
	eps=solver_tol

!	-Intrinsic GMRES constants
	gamma1 = -1
      
!-Set up preconditioner matrix (PD: See Heath's book.
!-Use Diagonal preconditioner)
	
	allocate(adiag(nsg),w(nsg),v(nsg,mmax),gamma(mmax))
	
	do i=1,nsg
	  adiag(i) = 1.0!-ain(i,i)
	enddo

	!Initialize
	call rzero(x,nsg)
	
	!Solve M*w=b
	call precon(w,b,adiag,nsg)
	
	!Compute L2 Norm
	gamma(1) = twonorm(w,nsg)
	
	if (gamma(1) == 0) return
	t = 1./gamma(1)
	call cmult2(v,w,t,nsg)

	if (gamma1 < 0) then
	  etol  = gamma(1)*eps
	  etol2 = gamma(1)*eps*eps
	else    
	  etol  = gamma1*eps     
	  etol2 = gamma1*eps*eps
	end if

      !Begin Arnoldi w/ modified gram-schmidt
	
	maxiter = mmax-1
	
	allocate(h(mmax,mmax),s(mmax),c(mmax))
	
      do k=1,maxiter
! 	write(*,*) k
         !construct v(:,k+1)=A*v(:,k)
         call lhs_gmres(v(:,k+1),v(:,k),ts,delta)

         !Solve M*w=v(:,k+1)
         call precon(w,v(:,k+1),adiag,nsg)

         do j=1,k
            h(j,k) = glsc2(w,v(:,j),nsg)
            t      = -h(j,k)
            call add2s2(w,v(:,j),t,nsg)
         end do !j

         !Compute L2 Norm
         h(k+1,k) = twonorm(w,nsg)

         if (abs(h(k+1,k)) <  etol2) then
            ifdone=.true.
         else
            ifdone=.false.
            t = 1./h(k+1,k)
            call cmult2(v(:,k+1),w,t,nsg)
         end if

         !apply Given's rotations to new column of H
         do i=1,k-1
            t = h(i,k)
            h(i  ,k) =  c(i)*t + s(i)*h(i+1,k)
            h(i+1,k) = -s(i)*t + c(i)*h(i+1,k)
         end do !i
         den        =  sqrt(  h(k,k)*h(k,k) + h(k+1,k)*h(k+1,k)  )
         c(k)       =  h(k  ,k) / den
         s(k)       =  h(k+1,k) / den
         h(k,k)     =  c(k)*h(k,k)+s(k)*h(k+1,k)
         gamma(k+1) = -s(k)*gamma(k)
         gamma(k  ) =  c(k)*gamma(k)
         if (ifdone .or. abs(gamma(k+1)) < etol) then
	      niter = k
	      exit
	   end if
      end do !k
	
	deallocate(adiag,w,s)
	
	!if we're here, we exceeded the max iteration, reduce k by 1
	if (k > maxiter) then
	  print*,' MAXITER exceeded'
	  k = k-1
	end if

      !Compute solution via back substitution
!     write(*,'("     k Tol     = ",i5,e16.8)')k,gamma(k+1)/gamma1
      kiter1=kiter0
      kiter0=k
      kiter=kiter + k
      do i=k,1,-1
         t = gamma(i)
         do j=k,i+1,-1
            t = t - h(i,j)*c(j)
         end do !j
         c(i) = t/h(i,i)
      end do !i

      !Sum up Arnoldi vectors
	do i=1,k
	  call add2s2(x,v(:,i),c(i),nsg)
	end do !i
	
	deallocate(v,h,c,gamma)
	
	end subroutine solve_gmres
!-----------------------------------------------------------------------
      subroutine add2s2(a,b,c1,n)
!     include 'param.h'
      real a(n), b(n)

      do i=1,n
        a(i)=a(i)+c1*b(i)
      end do  

      end subroutine add2s2
!-----------------------------------------------------------------------
      function twonorm(x,n)
!     include 'param.h'
      real x(n), twonorm

      tscal = x(1)*x(1)
      do i=2,n
         tscal = tscal+x(i)*x(i)
      end do
      twonorm=tscal
      if (twonorm.gt.0) twonorm = sqrt(twonorm)

      end function twonorm
!-----------------------------------------------------------------------
      function glsc2(x,y,n)
!     include 'param.h'
      real x(n), y(n), glsc2

      tscal = x(1)*y(1)
      do i=2,n
         tscal = tscal+x(i)*y(i)
      end do
      glsc2=tscal

      end function glsc2
!-----------------------------------------------------------------------
      subroutine cmult2(a,b,c,n)
!     include 'param.h'
      real a(n), b(n), c

      do i = 1, n
         a(i) = b(i)*c
      end do

      end subroutine cmult2
!-----------------------------------------------------------------------
      subroutine rzero(a,n)
!     include 'param.h'
      real a(n)

      do i = 1, n
         a(i) = 0.0
      end do   

      end subroutine rzero
!-----------------------------------------------------------------------
      subroutine precon(y,x,a,n)
!     include 'param.h'
      real y(n), x(n), a(n)

      do i=1,n
         y(i)=x(i)/a(i)
      end do !i

      end subroutine precon
!-----------------------------------------------------------------------
