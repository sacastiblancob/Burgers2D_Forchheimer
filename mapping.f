        subroutine mapping
! 	--------------------------------------------------------------
! 	This subroutine mapp the gllp into the physical space
! 	Jorge Escobar - Cornell University
! 	LAST REVISION 02/05/08
! 	--------------------------------------------------------------
! 	MODULES
		USE scrotum
		USE geom
		USE legendre
		USE map
		
	! 	LIST OF VARIABLES -> See modules
		
		implicit none

	! 	Local Variables
		integer :: c1,c2,ii,i,k,j
		real :: lx1,lx2,lz1,lz2

		do ii = 0,nsg-1,ns
		k = 0
		do i = 0,ns-1,n
			k = k+1
			do j = 1,n
			cx(ii+i+j) = points(j)
			cz(ii+i+j) = points(k)
			enddo
		enddo
		enddo

	! Los dos DO tenian nsubx y nsubz -1
		
		do k = 0,nsubz-1
		do j = 0,nsubx-1
			c1 = (k*nsubx) + (j+1)
			lx1 = (cgp(scp(c1,1),1))
			lx2 = (cgp(scp(c1,2),1))
			lz1 = (cgp(scp(c1,2),2))
			lz2 = (cgp(scp(c1,3),2))
			do i = 1,ns
			c2 = ns*((k*nsubx)+j) + i
			cx(c2) = ((lx2-lx1) * (cx(c2)+1)/2) + lx1
			cz(c2) = ((lz2-lz1) * (cz(c2)+1)/2) + lz1
			enddo
		enddo
		enddo
		
		end subroutine mapping
