	
	program Burgers
! 	-------------------------------------------------------
! 	This program solves the 2D Burgers' Equation 
! 	via Collocation (pseudo-spectral) 
! 	multidomain penalty method.

! 	Jorge Escobar - Cornell University
! 	Environmental Fluid Mechanics & Hydrology
! 	School of Civil and Environmental Engineering
! 	February - 2008
! 	--------------------------------------------------------
! 	MODULES
	USE scrotum
	USE mound
	USE geom
	USE legendre
	USE map
	USE aetas
	
	implicit none
	
	real,allocatable,dimension(:) :: ue,ve
	
	real, dimension(nsg) :: u,w
	real, allocatable, dimension(:) :: dudx,dwdx
	real, allocatable, dimension(:) :: dudz,dwdz
	real, allocatable, dimension(:) :: stx,stz
	real, allocatable, dimension(:) :: um1,um2,stxm1,stxm2
	real, allocatable, dimension(:) :: wm1,wm2,stzm1,stzm2
	real, allocatable, dimension(:) :: BGx,BGz
	real, dimension(n,n) :: F
	real, dimension(n) :: wg
	real, dimension(nsg) :: cont ! Vector verificar cont
	
	integer :: t, ff, ll, fn, pp, i, kk, ii, jj, temp11
	integer :: TotTimei,TotTimef
	real :: delta
	
	real, dimension(nsg) :: erru, errw 
	real, dimension(n*nsubx) :: bu, bw
	integer :: niterx, niterz
	real :: coux,couz, un2, wn2, filtro1
	real :: Linfu, Linfw, Lu, Lw, eu, ew
	
! 	Reading data file
	call readdata
	
! 	Generating Gauss-Lobatto-Legendre points
	allocate(points(n))
	call jacobl(pd, 0., 0.,points, n)

! 	Generating differentiation matrices
	allocate(d(n,n),d2(n,n),d3(n,n))
	call derv(pd,points,d,d2,d3,pd)

	filtro1 = 50.
	
! 	Generate Weights for Legendre polynomials and filter matrix
!	The value in line 60 (last of localfil arguments) is the filter value
	call quad(pd,points,wg,pd)
	call localFil(n, ns, points, F, wg, filtro1)
	
! 	Mapping from local coordinates to global coordinates
	allocate(cx(nsg),cz(nsg))
	call mapping

!	Esta parte imprime coordenadas y ya (APR) (190715)
!	open(33, file = "Coordinates.dat")
!	do i = 1,nsg
!            write(*,*) cx(i), cz(i)
!	enddo   
!	close(33)
	
!	stop
!	pause
	
! 	Setting exact solution
	allocate(ue(nsg),ve(nsg))
	
	call ex2dbur(ue, ve)
	
	fn = 10
	erru = 0.
	errw = 0.

!	call initburgers(ue,ve,u,w)
!	DESCOMENTAR LA ANTERIOR SI HAY PROBLEMAS LA DE ABAJO ES LA NUEVA
!	call initburgers(ue,ve, u, w)
	call ex2dbur(u, w)

	call output2(fn,u,w,erru,errw)

!	stop

! 	LLAMANDO EL PRIMER PASO DE VELOCIDAD PARA GENERAR CONDICIONES INICIALES
!	QUE NO SEAN CERO (SI SON CERO EL CFL ES 0.0) APR (170228)
	allocate(velocidades(n*nsubx,2))
	call readBvel(1)

! 	stop
! 	Setting Delta T and the maximum number of iterations
	call timestep(u,w)
! 	call CFL(u,w,coux,couz)
! 	write(*,*) coux,couz
! 	stop
	deallocate(velocidades)
	deallocate(points,d3)
	allocate(um1(nsg),um2(nsg),stxm1(nsg),stxm2(nsg))
	allocate(wm1(nsg),wm2(nsg),stzm1(nsg),stzm2(nsg))
	
	Lu = 10.
	Lw = 10.
	erru = 100.
	errw = 100.
!	pp = tmax

!	Aca se calcula el numero de pasos de tiempo y se imprime  (APR)
	pp = int(tmax/dT)
	write(*,*) 'Number of timesteps', pp
	write(*,*) 'Maximum time (seconds)', tmax
	write(*,*) 'Used dT for the simulation' , dT

!	El Sleep funciona para detener la corrida por el tiempo que se le 
!	diga entre parentesis (me ahorro poner go) (APR 170228)
!	call sleep(3)	
!	pause 

! 	Esta parte controla la cantidad de archivos que escribe el programa
!	Mover para obtener mas o menos refinamiento en los textos de resultado
	ll = pp / 1000
	ff = 1.0

!	Se localza el vector de velocidades (APR) - 170228
	allocate(velocidades(n * nsubx, 2))
	
	do t = 1,pp

!	  Llamando rutina para leer arreglo con top BC (170228)
	  call readBvel(t)
	 
! 	  write(*,*) t, pp
! 	  Setting the advective part

	  allocate(dudx(nsg),dwdx(nsg),dudz(nsg),dwdz(nsg))
	  
	  call diffx(u,dudx)
	  call diffx(w,dwdx)
	  
	  call diffz(u,dudz)
	  call diffz(w,dwdz)
	  allocate(stx(nsg),stz(nsg))
	  
! 	  Set advective term for X and Z-momentum equations
	  call spamer(u,w,dudx,dudz,dwdx,dwdz,stx,stz)

! 	  Setting the penalized BC - MODFIED APR (170228)
	  call BCPen(u,w,u,stx,ue,1) ! For X momentum eqn
	  call BCPen(u,w,w,stz,ve,2) ! For Z momentum eqn

!	ESTE PEDAZO ES EL ORIGINAL DEL PROGRAMA
!	Arriba quedan las modificaciones hechas para imponer BC	  
!! 	  Setting the penalized BC
!	  call BCPen(u,w,u,stx,ue) ! For X momentum eqn
!	  call BCPen(u,w,w,stz,ve) ! For Z momentum eqn
	  
! 	  Setting the penalized patching conditions
	  call patchpen(u,w,u,stx) ! For X momentum eqn
	  call patchpen(u,w,w,stz) ! For Z momentum eqn
	  
! 	  Advancing in time
	  call BDAB(t,u,stx,um1,um2,stxm1,stxm2)
	  call BDAB(t,w,stz,wm1,wm2,stzm1,stzm2)

!	  call output3(t,velocidades,bu,bw)
	
	  deallocate(stx,stz)

!	  Inclusion de filtros APR
! 	  Filtering Velocities
	  call filtering(n,numsub,ns,nsg,u,F)
	  call filtering(n,numsub,ns,nsg,w,F)
	  call interavg2d(t,u)
	  call interavg2d(t,w)

! 	  Estimating difference between BC and top velocities (as norm2)
!     Euclidean norm (190703 APR)	
	  jj = size(bu)
	  do kk = 0, nsubx - 1
		
		do ii = 1, n
		  temp11 = nsg - (ns * kk) + 1 - ii
		  !temp11 = 100 - (kk + 1) * n + ii
		  bu(jj) = u(temp11)
		  bw(jj) = w(temp11)		  
		  jj = jj - 1
		enddo
	  enddo
	
!	  stop 

	  un2 = norm2(velocidades(:,1) - bu)
	  wn2 = norm2(velocidades(:,2) - bw)

!	  call output3(t,velocidades,bu,bw)

!     Estimating norms of different vectors to check wether the BC are well or 
!	  ill imposed (190702 - APR)
!      write(*,*) "AFTER ADVECTIVE PART SOLUTION"
!	  write(*,*) "Norm of u difference: ", un2
!	  write(*,*) "Norm of w difference: ", wn2
!	  write(*,*) " "

!	  call sleep(1)
	 
! 	 Setting the implicit pressure treatment
	 
	  deallocate(dudx,dwdx,dudz,dwdz)
	   
! 	 Setting the diffusive part
	  
	  allocate(BGx(nsg),BGz(nsg))
	  
	  call setdelta(t,u,BGx,delta)
	  call setdelta(t,w,BGz,delta)

! 	  Imposing penalized BC in the RHS - APR (170228)
	  call BCrhs(delta,BGx,ue,1)
	  call BCrhs(delta,BGz,ve,2)
	
!! 	  Imposing penalized BC in the RHS
!	  call BCrhs(delta,BGx,ue)
!	  call BCrhs(delta,BGz,ve)
	  
! 	  Solving the system of equations
	  call solve_gmres(u,BGx,t,delta,niterx)
	  call solve_gmres(w,BGz,t,delta,niterz)

! 	  Estimating difference between BC and top velocities (as norm2)
!     Euclidean norm (190703 APR)
	  jj = size(bu)
	  do kk = 0, nsubx - 1		
		do ii = 1, n
		  temp11 = nsg - (ns * kk) + 1 - ii
		  !temp11 = 100 - (kk + 1) * n + ii
		  bu(jj) = u(temp11)
		  bw(jj) = w(temp11)		  
		  jj = jj - 1
		enddo
	  enddo

	  un2 = norm2(velocidades(:,1) - bu)
	  wn2 = norm2(velocidades(:,2) - bw)

!     Estimating norms of different vectors to check wether the BC are well or 
!	  ill imposed (190702 - APR)
      write(*,*) "AFTER DIFFUSIVE PART SOLUTION"
	  write(*,*) "Norm of u difference: ", un2
	  write(*,*) "Norm of w difference: ", wn2
	  write(*,*) " "

! 	  call sleep(1)

!  ! 	  Filtering Velocities
!	  call filtering(n,numsub,ns,nsg,u,F)
!	  call filtering(n,numsub,ns,nsg,w,F)
!	  call interavg2d(t,u)
!	  call interavg2d(t,w)
!     Printing differences between vectors	  
	  deallocate(BGx,BGz)

!     Verificando continuidad 
!	  cont = 0.
!	  call CFL(u,w,coux,couz)
!	  allocate(dudx(nsg),dwdz(nsg))
!	  call diffx(u,dudx)
!	  call diffz(w,dwdz)
!	  cont = dudx + dwdz
!	  write(*,*) 'CONTINUIDAD......', maxval(cont)

!	  deallocate(dudx,dwdz)
!	  call sleep(1)
	  
!	Comento lo que tenga que ver con los calculos de error porque no tienen 
!	nada que hacer en este caso (APR - 170228)
!	  call error(t,u,w,ue,ve,Linfu,Linfw,erru,errw)	  
!	  eu = abs(Linfu - Lu)/Linfu ! Stopping Criteria
!	  ew = abs(Linfw - Lw)/Linfu
	  
! 	  write(*,'(1X,I4,2F7.3,4E15.4)') t,coux,couz,eu,ew,Linfu,Linfw
!	  write(*,*) t, niterx, niterz
	  write(*,*) 'Tiempo REAL simulado:', t*dT
	  
!	  if (Linfu < 1e-5 .and. Linfw < 1e-5) then
!	   if (eu < 1e-8 .and. ew < 1e-8) then
!	    write(*,*) 'After ',t,'time steps'
!	    write(*,*) 'Linf norm for u = ', Linfu
!	    write(*,*) 'Linf norm for w = ', Linfw
!	    fn = fn + 1
!	    call system('cd test')
!	    call output2(fn,u,w)
!	    stop
!	   endif
!	  endif
	  
	  Lu = Linfu
	  Lw = Linfw
	    
	  if (t == ll * ff) then
	    ff = ff + 1.0
	    fn = fn + 1
		call output2(fn,u,w,erru,errw)
		call output3(fn,velocidades,bu,bw)
!		Escribir velocidad vertical!		
!		write(*,*) w
! 	    call error(t,u,w,ue,ve,Linfu,Linfw)
! 	    write(*,'(1X,I4,F7.3,3E15.4)') t, nu, dT*t, coux, couz
	  endif
	  
	enddo
	
	deallocate(d,d2)
	deallocate(cx,cz)
	deallocate(um1,um2,stxm1,stxm2)
	deallocate(wm1,wm2,stzm1,stzm2)
	deallocate(ue,ve)
	
!	Deslocalizando arreglo de velocidades
	deallocate(velocidades)
	
!	call timecall(TotTimef) Comentado APR
!	call timeelapse(TotTimei,TotTimef)  Comentado APR
	
	end program Burgers