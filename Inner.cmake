cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(shape4D)

# External dependencies
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH} )

find_package (FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIR})

#set (CMAKE_CXX_FLAGS "-fopenmp")

#find_package(MPI REQUIRED)

#SET(CMAKE_C_COMPILER mpicc)
#SET(CMAKE_CXX_COMPILER mpicxx)

#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

# Include directories
include_directories(
  include/
  src/
  #${MPI_INCLUDE_PATH}
)

FILE(GLOB_RECURSE shape4D_INCLUDE "include/*.h")

set(shape4D_SOURCE 	
   src/array1d.txx		
   src/array2d.txx
   src/array3d.txx
   src/adaptivegradientdescent.cpp
   src/algorithm.cpp
   src/shape4dstate.cpp
   src/grid.cpp
   src/gridoptimize.cpp
   src/helper.cpp
   src/landmarks.cpp
   src/optimizer.cpp
   src/regressionacceleration.cpp
   src/regression.cpp
   src/regressionparams.cpp
   src/regressionvelocity.cpp	
   src/runexperiment.cpp
   src/saveshapesandvectors.cpp
   src/surfacecurrent.cpp
   src/targetdata.cpp
   src/shapeobject.cpp
   src/tmplandmark.cpp
   src/tmpsurfacecurrent.cpp
   src/multiobjectcomplex.cpp
   src/vtkpolydatareader.cpp
   src/vtkpolydatawriter.cpp
   src/tinystr.cpp
   src/tinyxml.cpp
   src/tinyxmlerror.cpp
   src/tinyxmlparser.cpp
   src/main.cpp
)

add_executable(shape4D ${shape4D_SOURCE})
target_link_libraries(shape4D ${FFTW_LIBRARIES})
add_custom_target(include SOURCES ${shape4D_INCLUDE})

# If install subdirectory is not set, set it to default value
if(NOT INSTALL_RUNTIME_DEST)
  set(INSTALL_RUNTIME_DEST bin)
endif()

install(TARGETS shape4D RUNTIME DESTINATION ${INSTALL_RUNTIME_DEST} COMPONENT Executables)
