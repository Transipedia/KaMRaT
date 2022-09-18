add_library(dataStruct
			contig_elem.cpp
			feature_elem.cpp
			merge_knot.cpp
			scorer.cpp
)
target_link_libraries(dataStruct PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(dataStruct PRIVATE mlpack ${ARMADILLO_LIBRARIES})
target_include_directories(dataStruct PUBLIC
		"${PROJECT_SOURCE_DIR}/src/data_struct"
		"${PROJECT_SOURCE_DIR}/src/utils/"
)
