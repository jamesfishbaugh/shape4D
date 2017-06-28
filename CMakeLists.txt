cmake_minimum_required(VERSION 3.1)

project(shape4d)

if(shape4d_BUILD_SLICER_EXTENSION)
  #-----------------------------------------------------------------------------
  # Extension meta-information
  set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/shape4d")
  set(EXTENSION_CATEGORY "Shape Analysis")
  set(EXTENSION_CONTRIBUTORS "Francois Budin (Kitware Inc.), James Fishbaugh (NYU)")
  set(EXTENSION_DESCRIPTION "shape4D is an extension for estimating continuous evolution from a discrete collection of shapes, designed to produce realistic anatomical trajectories.")
  set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/shape4d.png")
  set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/shape4d/Screenshots/1.png")
  set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any
  set(EXTENSION_BUILD_SUBDIRECTORY inner-build)

  set(SUPERBUILD_TOPLEVEL_PROJECT inner)

  #-----------------------------------------------------------------------------
  # Extension dependencies
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
  mark_as_superbuild(Slicer_DIR)

  find_package(Git REQUIRED)
  mark_as_superbuild(GIT_EXECUTABLE)
  
  set(${EXTENSION_NAME}_SUPERBUILD ON CACHE BOOL "Build ${EXTENSION_NAME} and the projects it depends on.")
endif()

#-----------------------------------------------------------------------------
# SuperBuild setup
option(${EXTENSION_NAME}_SUPERBUILD "Build ${EXTENSION_NAME} and the projects it depends on." OFF)
mark_as_advanced(${EXTENSION_NAME}_SUPERBUILD)
if(${EXTENSION_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
# Shape4d application
include("${CMAKE_CURRENT_SOURCE_DIR}/Inner.cmake")

