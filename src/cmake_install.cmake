# Install script for directory: /Users/nlin/Documents/LEPH/ComputationalLib/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/ComputationalLib-0.2" TYPE STATIC_LIBRARY FILES "/Users/nlin/Documents/LEPH/ComputationalLib/src/libComputationalLib.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ComputationalLib-0.2/libComputationalLib.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ComputationalLib-0.2/libComputationalLib.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/ComputationalLib-0.2/libComputationalLib.a")
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ComputationalLib-0.2" TYPE FILE FILES
    "/Users/nlin/Documents/LEPH/ComputationalLib/include/ComputationalLib/csv_linreg.hpp"
    "/Users/nlin/Documents/LEPH/ComputationalLib/include/ComputationalLib/linreg.hpp"
    "/Users/nlin/Documents/LEPH/ComputationalLib/include/ComputationalLib/multiquad.hpp"
    "/Users/nlin/Documents/LEPH/ComputationalLib/include/ComputationalLib/poly.hpp"
    "/Users/nlin/Documents/LEPH/ComputationalLib/include/ComputationalLib/polyreg.hpp"
    "/Users/nlin/Documents/LEPH/ComputationalLib/include/ComputationalLib/qp.hpp"
    )
endif()

