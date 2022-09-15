cmake_minimum_required(VERSION 3.23)
project(DFT LANGUAGES CXX VERSION 0.1)

#-------------------------------------------------------------------------------
# Options

option( BUILD_SHARED_LIBS "Build shared libraries" true )

set(slate_DIR /opt/slate)
set(blaspp_DIR ${slate_DIR}/lib/blaspp)
set(lapackpp_DIR ${slate_DIR}/lib/lapackpp)

set( MPI_DETERMINE_LIBRARY_VERSION true )
find_package( MPI COMPONENTS CXX QUIET)
find_package( blaspp    REQUIRED QUIET)
find_package( lapackpp  REQUIRED QUIET)
find_package( slate     REQUIRED )

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set( prefix "/opt/DFT" CACHE PATH "Shorthand for CMAKE_INSTALL_PREFIX" )
    set( CMAKE_INSTALL_PREFIX "${prefix}"
            CACHE PATH
            "Install path prefix, prepended onto install directories."
            FORCE
            )
    message( STATUS "Setting CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
    list( APPEND CMAKE_SYSTEM_PREFIX_PATH ${CMAKE_INSTALL_PREFIX} )
else()
    message( STATUS "Using CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
endif()

# GNU Filesystem Conventions
include( GNUInstallDirs )

# the RPATH to be used when installing, but only if it's not a system directory
list( FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
        "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir )
if ("${isSystemDir}" STREQUAL "-1")
    list( APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" )
endif()

message(DEBUG "Settings:
CMAKE_INSTALL_PREFIX   = ${CMAKE_INSTALL_PREFIX}
CMAKE_BUILD_TYPE       = ${CMAKE_BUILD_TYPE}
BUILD_SHARED_LIBS      = ${BUILD_SHARED_LIBS}
slate_DIR              = ${slate_DIR}
blaspp_DIR             = ${blaspp_DIR}
lapackpp_DIR           = ${lapackpp_DIR}
" )

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
    cmake --build . #Runs the executable" )
endif()

#-------------------------------------------------------------------------------
# Build library.

file(
        GLOB_RECURSE libDFT_src
        CONFIGURE_DEPENDS
        src/*.cpp
)

message(DEBUG "libDFT_src = ${libDFT_src}" )
add_library(DFT ${libDFT_src})
target_link_libraries(DFT PUBLIC slate)

#-------------------------------------------------------------------------------
# Include directories
target_include_directories(
        DFT
        PUBLIC
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
            "$<INSTALL_INTERFACE:include>"
        PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

# Set compile features
set( gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU>" )
target_compile_features( DFT PUBLIC cxx_std_17)
target_compile_options(
        DFT PRIVATE "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall>>" )
set_target_properties(
        DFT PROPERTIES
        CXX_STANDARD_REQUIRED true  # prohibit < c++17
        CXX_EXTENSIONS false        # prohibit gnu++17
)

#-------------------------------------------------------------------------------
# Add 'make lib' target.
add_custom_target( lib DEPENDS DFT )

#-------------------------------------------------------------------------------
set(install_configdir "${CMAKE_INSTALL_LIBDIR}/DFT" )

# Install library and add to <package>Targets.cmake
install(
        TARGETS DFT
        EXPORT DFTTargets
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

# Install header files
install(
        # / copies contents, not directory itself
        DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        FILES_MATCHING REGEX "\\.(hpp|ipp)$"
)

# Install <package>Targets.cmake
install(
        EXPORT DFTTargets
        DESTINATION "${install_configdir}"
)

# Also export <package>Targets.cmake in build directory
export(
        EXPORT DFTTargets
        FILE "DFTTargets.cmake"
)

# Install <package>Config.cmake and <package>ConfigVersion.cmake,
# to enable find_package( <package> ).
include( CMakePackageConfigHelpers )
configure_package_config_file(
        "DFTConfig.cmake.in"
        "DFTConfig.cmake"
        INSTALL_DESTINATION "${install_configdir}"
)
write_basic_package_version_file(
        "DFTConfigVersion.cmake"
        VERSION "${DFT_VERSION}"
        COMPATIBILITY AnyNewerVersion
)
install(
        FILES "${CMAKE_CURRENT_BINARY_DIR}/DFTConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/DFTConfigVersion.cmake"
        DESTINATION "${install_configdir}"
)
