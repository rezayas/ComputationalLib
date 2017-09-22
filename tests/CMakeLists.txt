# Require C++14 compilation
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

include_directories (${TEST_SOURCE_DIR}/src
                     )

# Eyal's tests

add_executable (LinearTest
				lineartest.cpp)

add_executable (MinTest
				mintest.cpp)

add_executable (PolyTest
				polytest.cpp)

# Catch tests

add_executable (Test
				tests-main.cpp
				tests-CalibrateSin.cpp
				tests-CalibrateSinNDimensional.cpp)

# CalibrateSinRunner

add_executable (CalibrateSinRunner
				CalibrateSinRunner.cpp)

target_link_libraries (LinearTest
					   ComputationalLib
					   Eigen3::Eigen)

target_link_libraries (MinTest
					   ComputationalLib
					   Eigen3::Eigen)

target_link_libraries (PolyTest
					   ComputationalLib
					   Eigen3::Eigen)

target_link_libraries (Test
					   Catch
					   ComputationalLib
					   Eigen3::Eigen)

target_link_libraries (CalibrateSinRunner
					   ComputationalLib
					   Eigen3::Eigen)

add_test (NAME MyTest COMMAND Test WORKING_DIRECTORY <${CMAKE_CURRENT_SOURCE_DIR}>)
