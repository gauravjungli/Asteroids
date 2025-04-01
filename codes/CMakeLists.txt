# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(Landslides)
# Option to select which project to build
#option(BUILD_LANDSLIDES_1D "Build the landslides_1D project" OFF)
#option(BUILD_LANDSLIDES_2D "Build the landslides_2D project" ON)

# Include the appropriate subdirectory based on the selected option
if(BUILD_LANDSLIDES_1D)
    add_subdirectory(${CMAKE_SOURCE_DIR}/landslides_1D)
endif()

if(BUILD_LANDSLIDES_2D)
    add_subdirectory(${CMAKE_SOURCE_DIR}/landslides_2D)
endif()
