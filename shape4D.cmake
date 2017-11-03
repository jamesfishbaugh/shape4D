cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

#-----------------------------------------------------------------------------
# Slicer extension
#-----------------------------------------------------------------------------
if(EXTENSION_SUPERBUILD_BINARY_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
# External dependencies
#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH} )

option(USE_VTK "Use VTK instead of built-in vtk polydata reader/writer" OFF)
option(USE_SEM "Use SlicerExecutionModel as CLI" OFF)

find_package (FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIR})

# If we build a Slicer CLI module
# This needs to be before `find_package(VTK REQUIRED)`
if(USE_SEM)
  add_definitions("-DUSE_SEM")
  find_package(SlicerExecutionModel REQUIRED)
  include(${SlicerExecutionModel_USE_FILE})
endif()

# If we use VTK file reader instead of custom file reader
# This needs to be before `SEMMacroBuildCLI(...)` to include
# ${VTK_LIBRARIES}
if(USE_VTK)
    add_definitions("-DUSE_VTK")
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
endif()

#-----------------------------------------------------------------------------
# shape4D application
#-----------------------------------------------------------------------------

# Include directories
include_directories(
  include/
  src/
  )

FILE(GLOB_RECURSE ${PROJECT_NAME}_INCLUDE "include/*.h")

set(${PROJECT_NAME}_SOURCE
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
  src/polydatareader.cpp
  src/polydatawriter.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
  src/main.cpp
  )

if(USE_SEM)
  # Build a Slicer CLI
  add_definitions("-DUSE_SEM")
  SEMMacroBuildCLI(
    NAME ${PROJECT_NAME}
    ADDITIONAL_SRCS ${${PROJECT_NAME}_SOURCE}
    TARGET_LIBRARIES ${FFTW_LIBRARIES} ${VTK_LIBRARIES}
    INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    )
else()
  # Build an independent executable
  add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCE})
  target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARIES} ${VTK_LIBRARIES})
endif()

# Show header files in IDE
add_custom_target(include SOURCES ${${PROJECT_NAME}_INCLUDE})

#-----------------------------------------------------------------------------
# Tests
#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(testing)
endif()


#-----------------------------------------------------------------------------
# Install
#-----------------------------------------------------------------------------

#install(TARGETS shape4D RUNTIME DESTINATION ${INSTALL_RUNTIME_DEST} COMPONENT Executables)
# Slicer extension packaging
if(EXTENSION_SUPERBUILD_BINARY_DIR)
  install(FILES ${FFTW_INSTALL_LIBRARIES} DESTINATION ${Slicer_INSTALL_THIRDPARTY_LIB_DIR})
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
  include(${Slicer_EXTENSION_CPACK})
endif()
