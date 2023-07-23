CXX = g++
CXXFLAGS  =  -O3 -std=c++17      
#CXXFLAGS :=  -g -std=c++17

objects =  main.o Omega.o bc.o Characteristics.o \
            Sphere.o gravity.o IC.o IO.o March.o \
			pressure_shed.o solver.o topography.o TVD.o FS.o
objects_test =  main.o Omega.o bc.o Characteristics.o \
            test.o gravity.o IC.o IO.o March.o \
			pressure_shed.o solver.o topography.o TVD.o	FS.o	
all: gaurav,test

test: $(objects_test)
	$(CXX) $(CXXFLAGS)  -o test $(objects_test)

gaurav: $(objects)
		$(CXX) $(CXXFLAGS) -o gaurav $(objects)
$(objects): gauravlib.h
$(objects_test): gauravlib.h

clean : 
		rm gaurav test $(objects)