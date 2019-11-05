# Set fortran compiler here
FC=gfortran


%.o : %.f90; $(FC) -c $(FFLAGS) $< -o $@

OBJECTS= get-allowed-3ph-y.o

SC-FERMI: $(OBJECTS)
	  rm -f get-allowed-3ph-y
	  $(FC) -o get-allowed-3ph-y $(OBJECTS) 
	  rm -f get-allowed-3ph-y.o

clean:
	-rm -f *.o; touch *.f90
