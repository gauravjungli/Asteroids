cmake_minimum_required(VERSION 3.10) # Ensure compatibility
project(Landslides_2D VERSION 1.0)              # Name your CMake project

# Compilation Settings (Similar to Makefile variables)
set(CMAKE_CXX_STANDARD 17)              # Set C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)     # Enforce the standard
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3") # Optimization level 

# Source Files (Similar to Makefile objects)
set(SOURCE_FILES 
     
     src/bc.cpp
     src/Characteristics.cpp
     src/Sphere.cpp
    src/gravity.cpp
     src/IC.cpp
    src/IO.cpp
    src/March.cpp
     src/pressure_shed.cpp
    src/solver.cpp
     src/topography.cpp
   src/TVD.cpp 
    src/FS.cpp
)

# Header File Dependency
set(HEADER_FILE src/gauravlib.h)

# Executables (Targets in CMake)
add_executable(landslide  src/main.cpp ${SOURCE_FILES} ${HEADER_FILE}) # Main program
#add_executable(test  test.cpp  ${SOURCE_FILES} ${HEADER_FILE}) # Test program




