target_sources(DFT
        PRIVATE Potentials.cpp
        PUBLIC FILE_SET HEADERS
        BASE_DIRS ${PROJECT_SOURCE_DIR}
        FILES Potentials.h)