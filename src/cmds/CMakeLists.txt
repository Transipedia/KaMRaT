add_library(kamratIndex kamratIndex.cpp)
target_link_libraries(kamratIndex PRIVATE indexLoading seqCoding boost_iostreams)
target_include_directories(kamratIndex PRIVATE "${PROJECT_SOURCE_DIR}/src/runinfo_files/")

add_library(kamratMerge kamratMerge.cpp)
target_link_libraries(kamratMerge PRIVATE indexLoading seqCoding vectOp dataStruct boost_iostreams)
target_include_directories(kamratMerge PRIVATE "${PROJECT_SOURCE_DIR}/src/data_struct/")
target_include_directories(kamratMerge PRIVATE "${PROJECT_SOURCE_DIR}/src/runinfo_files/")

add_library(kamratFilter kamratFilter.cpp)
target_link_libraries(kamratFilter PRIVATE indexLoading)
target_include_directories(kamratFilter PRIVATE "${PROJECT_SOURCE_DIR}/src/runinfo_files/")

