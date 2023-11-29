add_executable(unittests
    unittests.cpp
    test_vect_opera.cpp
    test_scorer.cpp
)

target_link_libraries(unittests
    PRIVATE
    vectOp
    dataStruct
    armadillo
    indexLoading
)

target_include_directories(unittests
    PUBLIC
    "${PROJECT_SOURCE_DIR}/src/tests/"
)
