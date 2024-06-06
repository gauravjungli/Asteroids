cmake_minimum_required(VERSION 3.10) # Ensure compatibility
project(Asteroids VERSION 1.0)              # Name your CMake project

# Compilation Settings (Similar to Makefile variables)
set(CMAKE_CXX_STANDARD 17)              # Set C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)     # Enforce the standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3") # Optimization level 

# Source Files (Similar to Makefile objects)
set(SOURCE_FILES 
    main.cpp
    Omega.cpp
    bc.cpp
    Characteristics.cpp
    Sphere.cpp
    gravity.cpp
    IC.cpp
    IO.cpp
    March.cpp
    pressure_shed.cpp
    solver.cpp
    topography.cpp
    TVD.cpp 
    FS.cpp
)

# Header File Dependency
set(HEADER_FILE gauravlib.h)

# Executables (Targets in CMake)
add_executable(gaurav ${SOURCE_FILES} ${HEADER_FILE}) # Main program
add_executable(test   main.cpp test.cpp ${SOURCE_FILES} ${HEADER_FILE}) # Test program

# Custom Clean Target (Optional)
add_custom_target(clean
    COMMAND rm -f gaurav test ${SOURCE_FILES} 
    COMMENT "Cleaning up..."
)
