cmake_minimum_required(VERSION 3.19)

project(KaMRaT VERSION 2.0
                  DESCRIPTION "KaMRaT is a tool for finding substrings with interesting properties in large NGS datasets"
                  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# add_compile_options(-Ofast)


# Include local installed libraries

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_library(Armadillo::Armadillo INTERFACE IMPORTED)
set_target_properties(
  Armadillo::Armadillo PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${ARMADILLO_INCLUDE_DIR}"
  INTERFACE_LINK_LIBRARIES "${ARMADILLO_LIBRARIES}"
)

find_package(OpenMP)
include_directories(${OPENMP_INCLUDE_DIR})

add_subdirectory(external/mlpack)

# Recursive exploration of the project
add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(tests)
