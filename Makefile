# Set fortran compiler here
FC=gfortran


%.o : %.f90; $(FC) -c $(FFLAGS) $< -o $@

OBJECTS= get-allowed-3ph.o

SC-FERMI: $(OBJECTS)
	  rm -f get-allowed-3ph
	  $(FC) -o get-allowed-3ph $(OBJECTS) 
	  rm -f get-allowed-3ph.o

clean:
	-rm -f *.o; touch *.f90
