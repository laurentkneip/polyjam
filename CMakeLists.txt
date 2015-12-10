cmake_minimum_required(VERSION 2.4.6)
project(polyjam)

set (POLYJAM_VERSION_MAJOR 1)
set (POLYJAM_VERSION_MINOR 0)

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

add_definitions (-Wall -march=native -O3 -msse3 -mssse3 -fomit-frame-pointer -std=c++11 )
cmake_policy(SET CMP0005 OLD)

IF(NOT DEFINED MACAULAYCOMMAND)
  SET(MACAULAYCOMMAND '"M2"')
#  SET(MACAULAYCOMMAND '"/./Applications/Macaulay2-1.7/bin/M2"')
ENDIF()
ADD_DEFINITIONS(-DMACAULAYCOMMAND=${MACAULAYCOMMAND})

SET(TEMPPATH '"${PROJECT_SOURCE_DIR}/polyjam_M2_scripts/"')
ADD_DEFINITIONS(-DTEMPPATH=${TEMPPATH})
SET(SOLVERPATH '"${PROJECT_SOURCE_DIR}/polyjam_solvers/"')
ADD_DEFINITIONS(-DSOLVERPATH=${SOLVERPATH})

#set( OpenCV_LIBS
#     /usr/lib/libopencv_core.so
#     /usr/lib/libopencv_highgui.so )

include_directories(${PROJECT_SOURCE_DIR}/include)

#get modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules/")

# get eigen
find_package( Eigen REQUIRED)
include_directories( ${EIGEN_INCLUDE_DIRS} )
# get OpenCV
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

set( POLYJAM_SOURCE_FILES
  src/polyjam.cpp
  src/fields/R.cpp
  src/fields/Q.cpp
  src/fields/Sym.cpp
  src/fields/Zp.cpp
  src/core/Coefficient.cpp
  src/core/Monomial.cpp
  src/core/Term.cpp
  src/core/Poly.cpp
  src/core/PolyMatrix.cpp
  src/elimination/CMatrix.cpp
  src/elimination/methods.cpp
  src/elimination/GaussJordan.cpp
  src/exports/MacaulayTwo.cpp )

set( POLYJAM_HEADER_FILES
  include/polyjam/polyjam.hpp
  include/polyjam/fields/Field.hpp
  include/polyjam/fields/R.hpp
  include/polyjam/fields/Q.hpp
  include/polyjam/fields/Sym.hpp
  include/polyjam/fields/Zp.hpp
  include/polyjam/core/Coefficient.hpp
  include/polyjam/core/Monomial.hpp
  include/polyjam/core/Term.hpp
  include/polyjam/core/Poly.hpp
  include/polyjam/core/PolyMatrix.hpp
  include/polyjam/elimination/GaussJordan.hpp
  include/polyjam/elimination/CMatrix.hpp
  include/polyjam/elimination/methods.hpp
  include/polyjam/exports/MacaulayTwo.hpp )

add_library( polyjam SHARED ${POLYJAM_SOURCE_FILES} ${POLYJAM_HEADER_FILES} )
target_link_libraries( polyjam ${OpenCV_LIBS} )

add_Executable( sw6pt polyjam_problems/sw6pt.cpp )
target_link_libraries( sw6pt polyjam )
