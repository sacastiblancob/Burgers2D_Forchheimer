        subroutine output3(fn, velocidades, bu, bw)
    ! 	----------------------------------------------------
    ! 	This subroutine writes boundary velocities imposed 
    !   and estimated by the 2D Burgers model
    ! 	Antonio Preziosi-Ribero
    ! 	Universidad Nacional de Colombia 
    ! 	July 2019
    ! 	----------------------------------------------------
    ! 	MODULES
        USE scrotum
        USE map
        
        implicit none
    ! 	Dummy arguments
        real, dimension(nsubx*n), intent(in) :: bu,bw
        real, dimension(nsubx*n, 2), intent(in) :: velocidades
        integer, intent(in) :: fn
    ! 	Local Variables
        integer :: k, tot 
        real :: c1,c2,c3,c4
        character (LEN=80) :: fout
        character (LEN=9 ) :: fmt2
        character (LEN=30) :: fspec
        
        fmt2='(i4.4,a4)'
        fspec = 'Boundary_'
        fout = fspec
        
        write(unit=fout(10:23),fmt=fmt2) fn,'.dat'
    
        open(fn,file=fout)  
    
        write(fn,*) 'Boundary values for 2D Burgers'
        write(fn,*) 'COLUMNS = "Imp u" "Imp w" "Calc u" "Calc w"'
    
    ! Total of points to ease things (nothing special)
        tot = n * nsubx
    
        do k = 0, tot - 1
          
            c1 = velocidades(k + 1, 1)
            c2 = velocidades(k + 1, 2)
            c3 = bu(k + 1)
            c4 = bw(k + 1)

        write(fn,'(*(F14.8))') c1, c2, c3, c4 

        enddo
        close(fn)
        
        end subroutine output3  
