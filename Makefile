#Compiler flags.
FC = mpiifort
FFLAGS = -O2 -fopenmp

#Created objects and modules.
OBJECTS = ISO_Precisions.o CSR_Type.o Nullarbor.o Spinifex.o Example.o
MODS = iso_precisions.mod csr_type.mod nullarbor.mod spinifex.mod

%.o: %.f90
	$(FC) $(FFLAGS) $< -c 

Example: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o Example 

Spinifex.o: ISO_Precisions.o CSR_Type.o Nullarbor.o Spinifex.f90

Nullarbor.o: ISO_Precisions.o CSR_Type.o Nullarbor.f90

CSR_Type.o: ISO_Precisions.o CSR_Type.f90

ISO_Precisions.o: ISO_Precisions.f90

.PHONY: test clean veryclean

#Clean .o and .mod from src	
clean:
	rm *.mod *.o
