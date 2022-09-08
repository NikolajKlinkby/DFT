target_sources(DFT
        PRIVATE Space1D.cpp
        PUBLIC FILE_SET HEADERS
        BASE_DIRS ${PROJECT_SOURCE_DIR}
        FILES Space1D.h)
add_subdirectory("Potentials")
add_subdirectory("Observables")
add_subdirectory("SCF")
add_subdirectory("Orbitals")

