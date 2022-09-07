cmake_minimum_required(VERSION 3.23)
project(example LANGUAGES CXX)

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
    make" )
endif()
#-------------------------------------------------------------------------------

set(DFT_DIR "../DFT-install/lib/cmake/DFT")
find_package(DFT)

add_executable(example)

target_sources(example PRIVATE "example.cpp")
target_link_libraries(example PRIVATE DFT::DFT)

get_target_property(DFT_dll DFT::DFT IMPORTED_LOCATION_RELEASE)

add_custom_target(mpi_example COMMAND mpirun -np 4 ./example)

