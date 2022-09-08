target_sources(DFT
        PRIVATE SCF.cpp
        PUBLIC FILE_SET HEADERS
        BASE_DIRS ${PROJECT_SOURCE_DIR}
        FILES SCF.h)