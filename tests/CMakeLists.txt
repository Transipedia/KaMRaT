add_executable(unittests
    test_vect_opera.cpp
)

target_link_libraries(unittests
    PRIVATE
    vectOp
)

target_include_directories(unittests
    PUBLIC
    "${PROJECT_SOURCE_DIR}/src/tests/"
)