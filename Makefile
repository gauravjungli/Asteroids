CXX = g++
CXXFLAGS := -g #-Wall -O2

objects =  main.o Angular_mom.o bc.o Characteristics.o \
            cv_flux_source.o gravity.o IC.o IO.o march.o \
			pressure_shed.o solver.o topography.o TVD.o
all: gaurav


gaurav: $(objects)
		$(CXX) $(CXXFLAGS) -o gaurav $(objects)
$(objects): gauravlib.h parameters.h

clean : 
		rm gaurav $(objects)