	subroutine readdata
! 	-----------------------------------------
! 	This subroutine reads the geometry file
! 	Jorge Escobar - Cornell University
! 	New version with MODULES
! 	LAST REVISION 02/04/2008
! 	-----------------------------------------
! 	MODULES
	USE scrotum
	USE geom

	implicit none

! 	Local variables
	integer :: i,j

	open(55,file = 'PruebaGeom.dat',status='old')

! 	Reading coordinates of grid points
	do i = 1,ngp
	  read(55, *)(cgp(i, j), j = 1, 2)
	enddo

! 	Reading Grid points of each subdomain
	do i = 1,numsub
	  read(55, *)(scp(i, j), j = 1, 4)
	enddo

! 	Reading Type of boundary conditions for each
! 	side of the square domain
	read(55,*)(cond(j), j = 1, 4)

! 	Reading the magnitude of the imposed BC
	read(55,*)(val(j), j = 1, 4)
	
	end subroutine readdata
