add_executable(kamrat
    kamrat.cpp
)

target_link_libraries(kamrat
    PRIVATE
    kamratIndex
    kamratFilter
    kamratMerge
    kamratScore
    kamratMask
    kamratQuery
)

target_include_directories(kamrat
    PUBLIC
    "${PROJECT_SOURCE_DIR}/src/include/"
)