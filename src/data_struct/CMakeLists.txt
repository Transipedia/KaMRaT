add_library(dataStruct
			contig_elem.cpp
			feature_elem.cpp
			merge_knot.cpp
)
target_include_directories(dataStruct PUBLIC "${PROJECT_SOURCE_DIR}/src/data_struct")