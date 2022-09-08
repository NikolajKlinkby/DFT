target_sources(DFT
        PRIVATE Observables.cpp
        PUBLIC FILE_SET HEADERS
        BASE_DIRS ${PROJECT_SOURCE_DIR}
        FILES Observables.h)