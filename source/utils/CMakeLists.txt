target_sources(DFT
        PRIVATE utils.cpp
        PUBLIC FILE_SET HEADERS
        BASE_DIRS ${PROJECT_SOURCE_DIR}
        FILES utils.h)