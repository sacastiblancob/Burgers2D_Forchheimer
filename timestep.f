	subroutine timestep(u,w)
! 	--------------------------------------------------
! 	This subroutine calculates the size of the timestep
! 	and the number of time steps to be used in the
! 	calculation

! 	Jorge Escobar - Cornell University
! 	February 2008
! 	--------------------------------------------------
! 	Variables
!	u, w -> x and z velocity fields
! 	tmax -> Max number of time steps
! 	dT -> Size of the time step

!	CN -> Courant Number
!	DeltaX -> Minimum grid spacing (in local coordinates)

	USE scrotum
	USE map
	USE aetas
	USE geom

	implicit none

! 	Dummy arguments
	real, dimension(nsg), intent(in) :: u, w

! 	Local variables
	real :: DeltaX, Ucn, Wcn, L
	real :: mindis, vel
	integer :: i, p
	real,dimension(n * nsubx) :: vel1, vel2

! 	L = 1.0 * 2.0 * acos(-1.0) ! For the circular bump
! 	L = 10.0 ! For the horizontal one
! 	L = 0.04 ! For Dam-break problem
! 	L = 200*24*3!60!*60
	L = 1.0!0.02

	mindis = 100.
	do i = 0,nsubx-1
	 p = (i*ns) + 1
	 DeltaX = abs(cx(p+1) - cx(p))
	 mindis = min(DeltaX,mindis)
	enddo

	do i = 0,nsubz-1
	 p = (i*ns*nsubx) + 1
	 DeltaX = abs(cz(p+n) - cz(p))
	 mindis = min(DeltaX,mindis)
	enddo

!	Asignando los valores de u y w a las columnas del arreglo
!	APR (170228)
	do i=1,nsubx*n
	  vel1(i) = abs(velocidades(i,1))
	enddo

	do i=1,nsubx*n
	  vel2(i) = abs(velocidades(i,2))
	enddo

!	SE CAMBIA PARA QUE SEA CONSISTENTE
!	Ucn = maxval(vel1)
!	Wcn = maxval(vel2)

	Ucn = maxval(u)
	Wcn = maxval(w)

	write(*,*) 'Maximas velocidades encontradas: ', Ucn, Wcn
	call sleep(3)

! 	write(*,*) maxval(u),maxval(w),CN

	if (Ucn > Wcn) then
	 vel = Ucn
	else
	 vel = Wcn
	endif

!	dT = CN * mindis / vel

!	dT = 1e-2

!       Timestep from Northwestern results
        dT = 3e-2

	CN = dT * vel / mindis

!       Esto ya no funciona (marzo/2020) deprecated
!	tmax = L / dT

	write(*,*) 'Delta t: ', dT
	write(*,*) 'Mindis y vel', mindis, vel
	write(*,*) 'timesteps: ', int(tmax/dT)
	write(*,*) 'CFL CALCULADO: ', CN

	end subroutine timestep
