add_library(kamratIndex kamratIndex.cpp)
target_link_libraries(kamratIndex PRIVATE indexLoading seqCoding boost_iostreams)
target_include_directories(kamratIndex PRIVATE "${PROJECT_SOURCE_DIR}/src/runinfo_files/")

add_library(kamratMerge kamratMerge.cpp)
target_link_libraries(kamratMerge PRIVATE indexLoading seqCoding vectOp dataStruct boost_iostreams mlpack-interface armadillo)
target_include_directories(kamratMerge PRIVATE "${PROJECT_SOURCE_DIR}/src/runinfo_files/")

add_library(kamratRank kamratRank.cpp)
target_link_libraries(kamratRank PRIVATE indexLoading dataStruct)
target_include_directories(kamratRank PRIVATE "${PROJECT_SOURCE_DIR}/src/runinfo_files/")

add_library(kamratFilter kamratFilter.cpp)
target_link_libraries(kamratFilter PRIVATE indexLoading)
target_include_directories(kamratFilter PRIVATE "${PROJECT_SOURCE_DIR}/src/runinfo_files/")

add_library(kamratMask kamratMask.cpp)
target_link_libraries(kamratMask PRIVATE indexLoading seqCoding)
target_include_directories(kamratMask PRIVATE "${PROJECT_SOURCE_DIR}/src/runinfo_files/")

add_library(kamratQuery kamratQuery.cpp)
target_link_libraries(kamratQuery PRIVATE indexLoading seqCoding)
target_include_directories(kamratQuery PRIVATE "${PROJECT_SOURCE_DIR}/src/runinfo_files/")

