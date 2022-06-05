add_library(indexLoading index_loading.cpp)
target_include_directories(indexLoading PUBLIC "${PROJECT_SOURCE_DIR}/src/utils/")

add_library(seqCoding seq_coding.cpp)
target_include_directories(seqCoding PUBLIC "${PROJECT_SOURCE_DIR}/src/utils/")

add_library(vectOp vect_opera.cpp)
target_include_directories(vectOp PUBLIC "${PROJECT_SOURCE_DIR}/src/utils/")
