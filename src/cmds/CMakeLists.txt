add_library(kamratFilter
    kamratFilter.cpp
)

target_link_libraries(kamratFilter
    PRIVATE
    indexLoading
)

target_include_directories(kamratFilter
    PRIVATE
    "${PROJECT_SOURCE_DIR}/src/runinfo_files/"
)
