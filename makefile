
#  Output file name
   EXC = Toni.out
   EXC2 = Toni2.out
# 
# loading routines to compile. 
#  OPTIONS here.   
#

# Se remueve la rutina error.o de la lista de rutinas que utiliza el programa
# error.o se encontraba despues de patching.o y antes de interavg2d.o (APR)
# Agregado readBvel.o antes de diffx.o (de acuerdo con orden de main)

USE_INCL= scrotum.o mound.o geom.o legendre.o map.o aetas.o\
	  main.o readdata.o gll.o derv.o mapping.o\
	  ex2dbur.o output2.o velocity.o timestep.o readBvel.o\
	  diffx.o diffz.o spamer.o BCPen.o patchpen.o BDAB.o\
	  setdelta.o initburgers.o CFL.o\
	  BCrhs.o solve_gmres.o lhs_gmres.o allmixed.o\
	  BC.o patching.o interavg2d.o quad.o localFil.o filtering.o\
	  output3.o 
	  
USE_ROUT= scrotum.f mound.f geom.f legendre.f map.f aetas.f\
	  main.f readdata.f gll.f derv.f mapping.f ex2dbur.f output2.f\
	  velocity.f timestep.f readBvel.f diffx.f diffz.f spamer.f BCPen.f\
	  patchpen.f BDAB.f setdelta.f initburgers.f CFL.f BCrhs.f solve_gmres.f\
	  lhs_gmres.f allmixed.f BC.f patching.f interavg2d.f quad.f localFil.f\
	  filtering.f output3.f
	      

#FC = /opt/intel/bin/ifort
FC = gfortran

#FFLAGS = -r8 -align none
FFLAGS = -O3 -fdefault-real-8 -ffixed-line-length-132

#  Compile

$(EXC): $(USE_INCL)
		$(FC) $(FFLAGS) $(USE_INCL) -o $@

$(USE_INCL):	$(INCLUDES)

# Estas lineas que siguen son la "receta" para el make all que pide Eclipse
# prueba hecha el 21 de abril de 2017
all:
	$(FC) $(FFLAGS) $(USE_ROUT) 

.f.o:
		$(FC) $(FFLAGS) -c $*.f

clean:
		rm -f core $(EXC) $(USE_INCL) *.mod
		echo Clean Done
