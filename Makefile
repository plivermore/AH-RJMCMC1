#Choose your favourite Fortran90 compiler here:
#F90=ifort
F90=gfortran -O3 -ffast-math -ffree-line-length-0


source=Fortran_source

all :  AH

sort.o : $(source)/sort.F90
	$(F90) -c $(source)/sort.F90
SUBS_AH.o : $(source)/SUBS_AH.F90 sort.o
	$(F90) -c $(source)/SUBS_AH.F90  
AH: $(source)/AH_RJMCMC.F90 SUBS_AH.o sort.o
	$(F90) -o AH $(source)/AH_RJMCMC.F90  SUBS_AH.o sort.o



.PHONY: clean

clean:
	rm *.o *.mod 
