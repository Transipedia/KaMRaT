cmake_minimum_required(VERSION 3.19)

project(KaMRaT VERSION 2.0
                  DESCRIPTION "KaMRaT is a tool for finding substrings with interesting properties in large NGS datasets"
                  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)


find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})


add_subdirectory(src)
add_subdirectory(apps)
