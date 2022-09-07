cmake_minimum_required(VERSION 3.23)
project(DFT LANGUAGES CXX VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

if (MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

set(CMAKE_PREFIX_PATH /opt/slate)

set( MPI_DETERMINE_LIBRARY_VERSION true )
find_package( MPI COMPONENTS CXX )
find_package( blaspp    REQUIRED )
find_package( lapackpp  REQUIRED )
find_package( slate     REQUIRED )

add_library(DFT SHARED)
target_link_libraries(DFT slate)
target_include_directories(DFT PRIVATE "${PROJECT_SOURCE_DIR}")
add_subdirectory("utils")

# Setting structure
install(TARGETS DFT
        EXPORT DFTTargets
        #installs the public headers in the include directory with their relative path
        FILE_SET HEADERS
        # for executables and dll on Win
        RUNTIME DESTINATION bin
        # shared libraries
        LIBRARY DESTINATION lib
        # for static libraries
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include)

# Setting local variables
install(EXPORT DFTTargets
        NAMESPACE DFT::
        DESTINATION lib/cmake/DFT)

# Writing the config files
include(CMakePackageConfigHelpers)
write_basic_package_version_file("DFTConfigVersion.cmake"
        VERSION ${DFT_VERSION}
        COMPATIBILITY AnyNewerVersion)

install(FILES "DFTConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/DFTConfigVersion.cmake"
        DESTINATION lib/cmake/DFT)
