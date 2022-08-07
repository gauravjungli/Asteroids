CXX= g++

objects =  main.o Angular_mom.o bc.o Characteristics.o \
            cv_flux_source.o gravity.o IC.o IO.o march.o \
			pressure_shed.o solver.o topography.o TVD.o
gaurav : $(objects)
		$(CXX) -o gaurav $(objects)
     $(objects) : gauravlib.h