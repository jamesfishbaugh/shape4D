cmake_minimum_required(VERSION 3.1)

project(shape4D)

#-----------------------------------------------------------------------------
# Standalone vs Slicer extension option
#-----------------------------------------------------------------------------

# This option should be named after the project name, it corresponds to the
# option set to ON when the project is build by the Slicer Extension build
# system.

set(_default OFF)
set(_reason "${PROJECT_NAME}_BUILD_SLICER_EXTENSION is ON")
if(NOT DEFINED ${PROJECT_NAME}_BUILD_SLICER_EXTENSION AND DEFINED Slicer_DIR)
  set(_default ON)
  set(_reason "Slicer_DIR is SET")
endif()
option(${PROJECT_NAME}_BUILD_SLICER_EXTENSION "Build as a Slicer Extension" ${_default})

set(_msg "Checking if building as a Slicer extension")
message(STATUS ${_msg})
if(${PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  message(STATUS "${_msg} - yes (${_reason})")
else()
  message(STATUS "${_msg} - no (${PROJECT_NAME}_BUILD_SLICER_EXTENSION is OFF)")
endif()


#-----------------------------------------------------------------------------
# Extension meta-information
#-----------------------------------------------------------------------------

if(${PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/shape4D")
  set(EXTENSION_CATEGORY "Shape Analysis")
  set(EXTENSION_CONTRIBUTORS "Francois Budin (Kitware Inc.), James Fishbaugh (NYU)")
  set(EXTENSION_DESCRIPTION "shape4D is an extension for estimating continuous evolution from a discrete collection of shapes, designed to produce realistic anatomical trajectories.")
  set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/shape4D.png")
  set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/shape4D/Screenshots/1.png")
  # Required in Superbuild.cmake
  set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any
  set(SUPERBUILD_TOPLEVEL_PROJECT ${PROJECT_NAME})
  set(EXTENSION_BUILD_SUBDIRECTORY ${PROJECT_NAME}-build) # Required for packaging purposes
endif()


#-----------------------------------------------------------------------------
# Superbuild
#-----------------------------------------------------------------------------

# Set to OFF by default, to avoid changing behavior (backward compatibility).
set(_default "OFF")
if(${PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  set(_default "ON")
endif()
option(${PROJECT_NAME}_SUPERBUILD "Build ${PROJECT_NAME} and the projects it depends on." ${_default})
mark_as_advanced(${PROJECT_NAME}_SUPERBUILD)
if(${PROJECT_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
# Shape4D-build
#-----------------------------------------------------------------------------
include("${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.cmake")
