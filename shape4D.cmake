
#-----------------------------------------------------------------------------
# Slicer extension
#-----------------------------------------------------------------------------
if(${PROJECT_NAME}_BUILD_SLICER_EXTENSION)
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
  find_package(SlicerExecutionModel REQUIRED)
  include(${SlicerExecutionModel_USE_FILE})
endif()

# If we use VTK file reader instead of custom file reader
# This needs to be before `SEMMacroBuildCLI(...)` to include
# ${VTK_LIBRARIES}
if(USE_VTK)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
endif()

#-----------------------------------------------------------------------------
# shape4D application
#-----------------------------------------------------------------------------

set(CMAKE_INCLUDE_CURRENT_DIR ON)

configure_file(
  include/shape4Dconfig.h.in
  shape4Dconfig.h
  )

# Include directories
include_directories(
  include/
  src/
  )

file(GLOB_RECURSE ${PROJECT_NAME}_INCLUDE "include/*.h")

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
  ${${PROJECT_NAME}_INCLUDE}
  )

if(USE_SEM)
  # Build a Slicer CLI
  SEMMacroBuildCLI(
    NAME ${PROJECT_NAME}
    ADDITIONAL_SRCS ${${PROJECT_NAME}_SOURCE}
    TARGET_LIBRARIES ${FFTW_LIBRARIES} ${VTK_LIBRARIES}
    )
  if(WIN32)
    add_custom_command(TARGET shape4D POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory ${Slicer_CLIMODULES_BIN_DIR}/$<CONFIG>
      COMMAND ${CMAKE_COMMAND} -E copy ${FFTW_INSTALL_LIBRARIES} ${Slicer_CLIMODULES_BIN_DIR}/$<CONFIG>
      )
  endif()
else()
  # Build an independent executable
  add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCE})
  target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARIES} ${VTK_LIBRARIES})
  if(WIN32)
    add_custom_command(TARGET shape4D POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/$<CONFIG>
      COMMAND ${CMAKE_COMMAND} -E copy ${FFTW_INCLUDE_DIR}/libfftw3-3.dll ${CMAKE_BINARY_DIR}/$<CONFIG>
      )
  endif()
endif()

#-----------------------------------------------------------------------------
# Tests
#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(testing)
endif()

#-----------------------------------------------------------------------------
# Slicer extension packaging
#-----------------------------------------------------------------------------
if(${PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  if(NOT APPLE)
    install(FILES ${FFTW_INSTALL_LIBRARIES} DESTINATION ${Slicer_INSTALL_THIRDPARTY_LIB_DIR} COMPONENT RuntimeLibraries)
  else()
    get_filename_component(FFTW_RUNTIME_DIRECTORY ${FFTW_LIB} DIRECTORY)
    set(${EXTENSION_NAME}_FIXUP_BUNDLE_LIBRARY_DIRECTORIES ${FFTW_RUNTIME_DIRECTORY} CACHE STRING "List of fixup bundle library directories" FORCE)

    set(${EXTENSION_NAME}_CUSTOM_CONFIG "
      set(${EXTENSION_NAME}_FIXUP_BUNDLE_LIBRARY_DIRECTORIES \"${FFTW_RUNTIME_DIRECTORY}\")
")
  endif()

  #-----------------------------------------------------------------------------
  set(EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS)
  set(${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS "${EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS}" CACHE STRING "List of external projects to install" FORCE)

  #-----------------------------------------------------------------------------
  list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
  list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS}")
  include(${Slicer_EXTENSION_GENERATE_CONFIG})
  include(${Slicer_EXTENSION_CPACK})
endif()
