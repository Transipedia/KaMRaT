cmake_minimum_required(VERSION 3.19)

project(KaMRaT VERSION 2.0
                  DESCRIPTION "KaMRaT is a tool for finding substrings with interesting properties in large NGS datasets"
                  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)


# Include local installed libraries

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})


# Recursive exploration of the project
add_subdirectory(src)
add_subdirectory(apps)
