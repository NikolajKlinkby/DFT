cmake_minimum_required(VERSION 3.23)
project(example LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

#-------------------------------------------------------------------------------
# Enforce out-of-source build
string( TOLOWER "${CMAKE_CURRENT_SOURCE_DIR}" source_dir )
string( TOLOWER "${CMAKE_CURRENT_BINARY_DIR}" binary_dir )
if ("${source_dir}" STREQUAL "${binary_dir}")
    message( FATAL_ERROR
            "Compiling with CMake requires an out-of-source build. To proceed:
    rm -rf CMakeCache.txt CMakeFiles/   # delete files in ${CMAKE_CURRENT_SOURCE_DIR}
    mkdir build
    cd build
    cmake ..
    cmake --build .
    cmake --install . --prefix ../../DFT-install" )
endif()
#-------------------------------------------------------------------------------

set(DFT_DIR "../DFT-install/lib/DFT")
find_package(DFT)

set( MPI_DETERMINE_LIBRARY_VERSION true )
find_package( MPI COMPONENTS CXX QUIET)

add_executable(example)

target_sources(example PRIVATE example.cpp)
target_link_libraries(example PUBLIC DFT)

# Run the executable
add_custom_target(mpi_example ALL COMMAND mpirun -np 4 ./example )
