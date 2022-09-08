target_sources(DFT
        PRIVATE Initial.cpp
        PUBLIC FILE_SET HEADERS
        BASE_DIRS ${PROJECT_SOURCE_DIR}
        FILES Initial.h)