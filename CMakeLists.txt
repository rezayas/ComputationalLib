cmake_minimum_required(VERSION 3.3)
project (ComputationalLib)

# Require C++14 compilation
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Define library version
set(ComputationalLib_VERSION_MAJOR 0 CACHE STRING "major version" FORCE)
set(ComputationalLib_VERSION_MINOR 2 CACHE STRING "minor version" FORCE)
set(ComputationalLib_VERSION ${ComputationalLib_VERSION_MAJOR}.${ComputationalLib_VERSION_MINOR} CACHE STRING "version" FORCE)

# Define install destinations
set(include_dest "include/ComputationalLib-${ComputationalLib_VERSION}")
set(main_lib_dest "lib/ComputationalLib-${ComputationalLib_VERSION}")
set(lib_dest ${main_lib_dest}/${CMAKE_BUILD_TYPE})

# Define options here

# Add subdirectories
add_subdirectory(src)

# # Enable testing
enable_testing()
add_subdirectory(tests)

install(FILES ComputationalLib-config.cmake DESTINATION ${main_lib_dest})
install(EXPORT ComputationalLib DESTINATION ${lib_dest})
