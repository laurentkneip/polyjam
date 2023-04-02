cmake_minimum_required(VERSION 2.4.6)
project(test_sw6pt)

# Set the build type.  Options are:
#
#  None (CMAKE_C_FLAGS or CMAKE_CXX_FLAGS used)
#  Debug (CMAKE_C_FLAGS_DEBUG or CMAKE_CXX_FLAGS_DEBUG)
#  Release (CMAKE_C_FLAGS_RELEASE or CMAKE_CXX_FLAGS_RELEASE)
#  RelWithDebInfo (CMAKE_C_FLAGS_RELWITHDEBINFO or CMAKE_CXX_FLAGS_RELWITHDEBINFO
#  MinSizeRel (CMAKE_C_FLAGS_MINSIZEREL or CMAKE_CXX_FLAGS_MINSIZEREL) 

set(CMAKE_BUILD_TYPE Release)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_definitions(-Wall -O3)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules/")
find_package(Eigen REQUIRED)

include_directories(${EIGEN_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIR}/unsupported ${PROJECT_SOURCE_DIR}/include)

add_executable( test_sw6pt test.cpp sw6pt.cpp )
