cmake_minimum_required(VERSION 3.1)

project(shape4d)

# Superbuild is automatically set to ON if shape4d is built as a Slicer extension.
# Otherwise, it is set to OFF by default, to avoid changing behavior (backward compatibility).
if(shape4d_BUILD_SLICER_EXTENSION)
    set(${PROJECT_NAME}_SUPERBUILD ON CACHE BOOL "Build ${PROJECT_NAME} and the projects it depends on.")
endif()

# SuperBuild setup (default to OFF, see explanation above)
option(${PROJECT_NAME}_SUPERBUILD "Build ${PROJECT_NAME} and the projects it depends on." OFF)
mark_as_advanced(${PROJECT_NAME}_SUPERBUILD)
if(${PROJECT_NAME}_SUPERBUILD)
  #-----------------------------------------------------------------------------
  # Extension meta-information
  set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/shape4d")
  set(EXTENSION_CATEGORY "Shape Analysis")
  set(EXTENSION_CONTRIBUTORS "Francois Budin (Kitware Inc.), James Fishbaugh (NYU)")
  set(EXTENSION_DESCRIPTION "shape4D is an extension for estimating continuous evolution from a discrete collection of shapes, designed to produce realistic anatomical trajectories.")
  set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/shape4d.png")
  set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/shape4d/Screenshots/1.png")
  set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any
  # Required in Superbuild.cmake
  set(EXTENSION_BUILD_SUBDIRECTORY inner-build)
  set(SUPERBUILD_TOPLEVEL_PROJECT inner)
  #-----------------------------------------------------------------------------
  # Extension dependencies
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
  mark_as_superbuild(Slicer_DIR)

  find_package(Git REQUIRED)
  mark_as_superbuild(GIT_EXECUTABLE)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
# Shape4d application
include("${CMAKE_CURRENT_SOURCE_DIR}/Inner.cmake")

