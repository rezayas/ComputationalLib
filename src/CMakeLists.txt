# Set headers
set(header_path "${ComputationalLib_SOURCE_DIR}/include/ComputationalLib")
set(header ${header_path}/csv_linreg.hpp
           ${header_path}/linreg.hpp
           ${header_path}/multiquad.hpp
           ${header_path}/poly.hpp
           ${header_path}/polyreg.hpp
           ${header_path}/qp.hpp
           )

# Require C++14 compilation
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Boost REQUIRED)
find_package(StatisticalDistributionsLib REQUIRED)
find_package(SimulationLib REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Set source files
set(src csv_linreg.cpp
        linreg.cpp
        multiquad.cpp
        poly.cpp
        polyreg.cpp
        qp.cpp
        CalibrateSin.cpp)

# Define library target
add_library(ComputationalLib ${header} ${src})

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${header_path})

target_link_libraries(ComputationalLib PUBLIC 
                      StatisticalDistributionsLib 
                      SimulationLib 
                      Eigen3::Eigen)

target_include_directories(ComputationalLib PUBLIC
  $<BUILD_INTERFACE:${ComputationalLib_SOURCE_DIR}/include> # for headers when building
  # $<BUILD_INTERFACE:${StatisticalDistributionsLib_INCLUDE_DIRS}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # for config_impl.hpp when building
  $<INSTALL_INTERFACE:${include_dest}> # for client in install mode
  $<INSTALL_INTERFACE:${lib_dest}>) # for config_impl.hpp in install mode

install(TARGETS ComputationalLib EXPORT ComputationalLib DESTINATION "${lib_dest}")
install(FILES ${header} DESTINATION "${include_dest}")
# install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config_impl.hpp DESTINATION "${include_dest}")
