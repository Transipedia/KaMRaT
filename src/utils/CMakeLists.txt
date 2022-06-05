add_library(indexLoading
    index_loading.cpp
)

target_include_directories(indexLoading
    PUBLIC
    "${PROJECT_SOURCE_DIR}/src/utils/"
)
