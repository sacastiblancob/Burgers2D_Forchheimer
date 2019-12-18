	subroutine readBvel(t)
    
	USE aetas
	USE geom


	implicit none
	
	integer, intent(in) :: t
	!	Cuando la variable viene de un modulo no se declara aca (APR)
!	real, dimension(n*nsubx,2), intent(inout) :: velocidades	

	integer :: i, j
	!character (LEN=12) :: froot = "Norm_vel.csv"
	!character (LEN=11) :: fpath = "Normalized/"
	character (LEN=12) :: froot = "Filt_vel.csv"
	character (LEN=9) :: fpath = "Filtered/"
	character (LEN=5) :: fnum_c
	!character (LEN=29) :: ftot ! Para nornmalizados
	character (LEN=27) :: ftot
	character (LEN=9) :: fmat = '(I5.5,A5)'
	

!	Convierto la variable a string para que quede lista
!	write(unit=fnum_c(1:5), fmt = fmat), t
	write(unit=fnum_c(1:5), fmt = fmat) t

	ftot = fpath // fnum_c // froot

!	Probando si el nombre del archivo esta organizado
!	write(*,*) ftot

!	Abro el archivo a leer para sacar los datos necesarios
	open(23, file=ftot, status='old')

!	Iteracion en arreglo bidimensional leyendo velocidades
	do i = 1, 100 

		read(23,*)(velocidades(i,j), j = 1, 2)		

!	Este pedazo escribe la velocidad leida para commprobar funcion
!		write(*,*) i, velocidades(i, 1), velocidades(i, 2) 
	enddo

	close(23)
	   
	end subroutine readBvel
