	subroutine Lectura(mediavel)


	USE aetas
	USE geom


	implicit none

	integer :: i, fnum
	real,intent(out) :: mediavel
	character (LEN=12) :: froot = "Filt_vel.csv"
	character (LEN=17) :: fname
	character (LEN=78) :: fpath = "/home/toni/Documents/2D_Comparison &
                    &/170227_2DBurgers/Python_scripts/Interp_Vel/"
	character (LEN=95) :: ftot
	
	

	open(55, file= 'Limpio_10.dat', status='old')
	write(*,*) pasos
	do i = 1, pasos 
		read(55,*) crudas(i)
!		write(*,*) 'Velocidad en ', 'i', velocidades(i)	[prueba en pan
!		talla]
	enddo
	close(55)
	
	mediavel = sum(crudas) / pasos
		   
	end subroutine Lectura
