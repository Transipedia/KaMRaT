# Define a library interface that holds mlpack stuff
add_library(mlpack-interface INTERFACE)

find_package(ENSMALLEN "2.10" QUIET) # Try to find ensmallen on system
if (TARGET ensmallen::ensmallen)     # Check for the imported target
  # Link the target to 'mlpack-interface'
  target_link_libraries(mlpack-interface INTERFACE ensmallen::ensmallen)
else() # If not found, fetch and build from github release

  # ensmallen and mlpack are not suited to be used as subprojects, i.e. using 'add_subdirectory'
  # on a submodule or using FetchContent. Ex:
  #
  # include(FetchContent)
  # FetchContent_Declare(
  #   ENSMALLEN
  #   GIT_REPOSITORY https://github.com/mlpack/ensmallen
  #   GIT_TAG 2.21.0
  #   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/mlpack_install
  # )
  # FetchContent_MakeAvailable(ENSMALLEN)
  # add_executable(main main.cpp)
  # target_link_libraries(main ensmallen::ensmallen)
  #
  # So here we use ExternalProject to build ensmallen/mlpack as separate projects
  # and install them at ${CMAKE_BINARY_DIR}/mlpack_install.

  include(ExternalProject)
  ExternalProject_Add(
    ENSMALLEN
    GIT_REPOSITORY https://github.com/mlpack/ensmallen
    GIT_TAG 2.21.0
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/mlpack_install
  )

  # Since ExternalProject, unlike FetchContent, does not provide cmake targets,
  # we manually configure the interface library using installed files
  # at ${CMAKE_BINARY_DIR}/mlpack_install.

  target_include_directories(mlpack-interface INTERFACE
    ${CMAKE_BINARY_DIR}/mlpack_install/include
  )

  # Ensure that ensmallen is ready before building something that depends on mlpack-interface
  add_dependencies(mlpack-interface ENSMALLEN)
endif()

# We use the same strategy for mlpack
find_package(MLPACK "4" QUIET)
if (TARGET mlpack::mlpack)
  target_link_libraries(mlpack-interface INTERFACE mlpack::mlpack)
else()
  include(ExternalProject)
  ExternalProject_Add(
    MLPACK
    GIT_REPOSITORY https://github.com/mlpack/mlpack
    GIT_TAG 4.3.0
    CMAKE_ARGS -DDOWNLOAD_DEPENDENCIES=ON -DBUILD_CLI_EXECUTABLES=OFF -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/mlpack_install
  )
  add_dependencies(mlpack-interface MLPACK)
  target_include_directories(mlpack-interface INTERFACE
    ${CMAKE_BINARY_DIR}/mlpack_install/include
  )
endif()

# Now mlpack can be linked using 'target_link_libraries(mytarget mlpack-interface)
