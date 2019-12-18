	subroutine BCPen(u,w,rho,st,bv,idc)
! 	-----------------------------------------------
! 	In this subroutine Penalized Dirichlet 
! 	conditions are imposed for Advection equation.
! 
! 	Jorge Escobar-Vargas
! 	Cornell University
! 	February 2008
! 	-----------------------------------------------

!	SE AGREGA LA VARIABLE idc QUE ES EL DENTIFICADOR DE LA COLUMNA
!	DE VELOCIDADES LEIDA PARA CADA UNO DE LOS CASOS. APR (170228)

! 	MODULES
	USE scrotum
	USE geom
	USE mound

	implicit none
	
! 	Dummy Variables
	real, dimension(nsg),intent(in) :: u,w,rho
	real, dimension(nsg),intent(in) :: bv ! Boundary Values
	real, dimension(nsg),intent(inout) :: st
	integer, intent(in) :: idc
	
! 	Local Variables
	integer :: k, i, temp, t1, temp1, SS
	real :: lz, lx, fac, temp2, temp3
	real :: tau, omega, alpha
	real :: mfb, mft, mfr, mfl

!	Poniendo el vaor de temp3 de acuerdo a la coordenada que se maneja 
!	(APR)
	if (idc == 1) then
		temp3 = val(1)
	else if (idc == 2) then
		temp3 = val(1)
	endif
	
!	Verificando si velocidades entro a la subrutina - APR (170228)
!	write(*,*) velocidades(1,1), velocidades(1,2)
!	call sleep(1)
	
! 	Square domain (Global)
	fac = 2.0  ! 2.0 Original value on the right (APR)
	omega = 2. / (pd * (pd + 1.))

! 	Conditions for the bottom
	   do k = 0, nsubx - 1
	     lz = abs(cgp(scp(k + 1, 2),2) - cgp(scp(k + 1, 3),2))
	     mfb = 2.0 / lz
	     do i = 1, n
	       temp = (ns * k) + i
	       if (w(temp) > 0.0) then

		alpha = w(temp)

		tau = fac * mfb / (2.0 * omega)

! La condicion de contorno es un solo valor si esta activo la primera linea
! La siguiente linea se va a cero el termino alpha*val(1) APR (170228)

 		st(temp) = st(temp) + (tau * (alpha * rho(temp) - alpha * temp3))
! 		st(temp) = st(temp) + tau*alpha*rho(temp)
!		st(temp) = st(temp) + (tau*(alpha*rho(temp) - alpha*bv(temp)))
	       endif
	     enddo
	   enddo

!	Variable SS that tracks size - counter for imposing BC
	   SS = size(velocidades(:,1))

!	Conditions for the top strongly imposed. (APR 190725)
!	This is for testing the boundary conditions with no penalization
	!    do k = 0, nsubx - 1
	! 	do i = 1, n
	! 	  temp = nsg - (ns * k) + 1 - i
	! 	  st(temp) = velocidades(SS, idc)
	! 	  SS = SS - 1
	! 	enddo
	!    enddo

! 	Conditions for the top
	   do k = 0, nsubx - 1
	     lz = abs(cgp(scp(numsub-k,2),2)-cgp(scp(numsub-k,3),2))
	     mft = 2.0 / lz
		 do i = 1, n
			
	       temp = nsg - (ns * k) + 1 - i
		   temp1 = SS
!		   temp1 = SS - (k + 1) * n + i

!	Testing the order of the imposed BC (190722 - APR)
!		   write(*,*) 'Esto es!!!'
!		   write(*,*) temp, temp1 

	       if (w(temp) < 0.0) then
			  alpha = abs(w(temp))
			  tau = (fac * mft) / (2.0 * omega)

!       La condicion de contorno es un solo valor si esta activa la 
!       primera linea
!!!		st(temp) = st(temp) + (tau*(alpha*rho(temp) - alpha*val(3)))

!  VERIFICAR SI ESTA LÍNEA QUE SIGUE DEBE SER INCLUIDA
!  AHORA SON DOS LÍNEAS (LAS QUE EMPIEZAN CON st(temp)... )
!  ACTIVAR Y DESACTIVAR.
!			  st(temp) = st(temp) + tau*alpha*rho(temp)  
!			  st(temp) = st(temp) - alpha*velocidades(temp1,idc)
		      temp2 = tau*alpha*(rho(temp)-velocidades(temp1,idc))
			  st(temp)=st(temp) + temp2
			
		   endif
		   
!		    temp1 = temp1 - 1
            SS = SS - 1

	     enddo
	   enddo

! 	Conditions for the right
	   do k = 1,nsubz
	     lx = abs(cgp(scp(k*nsubx,2),1) - cgp(scp(k*nsubx,1),1))
	     mfr = 2.0 / lx
	     do i = 1,n
	       temp = ns*((k*nsubx)-1)+(i*n)
	       if (u(temp) < 0.0) then
		alpha = abs(u(temp))
		tau = fac * mfr / (2.0*omega)
!               Lo mismo que en las BC anteriores
! 		st(temp) = st(temp) + (tau*(alpha*rho(temp) - alpha*val(2)))
!		st(temp) = st(temp) + (tau*(alpha*rho(temp) - alpha*bv(temp)))
!       Lineas traidas del codigo de Jorge Escbar
		t1 = ns*((k-1)*nsubx)+((i-1)*n)+1 !Left
		st(temp) = st(temp) + (tau*(alpha*rho(temp) - alpha*rho(t1)))
	       endif
	     enddo
	   enddo

! 	Conditions for the left
	   do k = 0,nsubz-1
!            Quite el valor abolsuto de lx y lo calcule en la sguiente linea APR
        lx = abs(cgp(scp((k*nsubx)+1,2),1) - cgp(scp((k*nsubx)+1,1),1))
!	     lx = abs(lx)
	     mfl = 2.0 / lx
	     do i = 0,n-1
	       temp = ns*(k*nsubx)+(i*n)+1
	       if (u(temp) > 0.0) then
		alpha = u(temp)
		tau = fac * mfl / (2.0*omega)
!       Lineas traidas del codigo de Jorge Escobar para Navier-Stokes
!               st(temp) = st(temp) + (tau*(alpha*rho(temp) - alpha*val(4)))
!		st(temp) = st(temp) + (tau*(alpha*rho(temp) - alpha*bv(temp)))
		t1 = ns*(((k+1)*nsubx)-1)+((i+1)*n) !Right
 		st(temp) = st(temp) + (tau*(alpha*rho(temp) - alpha*rho(t1)))

	       endif
	     enddo
	   enddo
	
	end subroutine BCPen
