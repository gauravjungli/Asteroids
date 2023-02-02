CXX = g++
CXXFLAGS := -g -std=c++17 #-Wall -O2

objects =  main.o Angular_mom.o bc.o Characteristics.o \
            cv_flux_source.o gravity.o IC.o IO.o march.o \
			pressure_shed.o solver.o topography.o TVD.o FS.o
objects_test =  main.o Angular_mom.o bc.o Characteristics.o \
            test.o gravity.o IC.o IO.o march.o \
			pressure_shed.o solver.o topography.o TVD.o	FS.o	
all: gaurav,test

test: $(objects_test)
	$(CXX) $(CXXFLAGS) -o test $(objects_test)

gaurav: $(objects)
		$(CXX) $(CXXFLAGS) -o gaurav $(objects)
$(objects): gauravlib.h
$(objects_test): gauravlib.h

clean : 
		rm gaurav $(objects)