# Require C++14 compilation
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories (${TEST_SOURCE_DIR}/src
                     )

add_executable (LinearTest
				lineartest.cpp)

add_executable (MinTest
				mintest.cpp)

add_executable (PolyTest
				polytest.cpp)

target_link_libraries (LinearTest
					   ComputationalLib
					   Eigen3::Eigen)

target_link_libraries (MinTest
					   ComputationalLib
					   Eigen3::Eigen)

target_link_libraries (PolyTest
					   ComputationalLib
					   Eigen3::Eigen)

add_test (NAME MyTest COMMAND Test WORKING_DIRECTORY <${CMAKE_CURRENT_SOURCE_DIR}>)
