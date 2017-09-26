cmake_minimum_required(VERSION 3.1)

project(shape4D)

#-----------------------------------------------------------------------------
# Slicer extension
#-----------------------------------------------------------------------------

# Superbuild is automatically set to ON in the outer-build if shape4D is built
# as a Slicer extension.
if(${PROJECT_NAME}_BUILD_SLICER_EXTENSION AND NOT EXTENSION_SUPERBUILD_BINARY_DIR)
  set(${PROJECT_NAME}_SUPERBUILD ON CACHE BOOL "Build ${PROJECT_NAME} and the projects it depends on.")
endif()

if(${PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  #-----------------------------------------------------------------------------
  # Extension meta-information
  set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/shape4D")
  set(EXTENSION_CATEGORY "Shape Analysis")
  set(EXTENSION_CONTRIBUTORS "Francois Budin (Kitware Inc.), James Fishbaugh (NYU)")
  set(EXTENSION_DESCRIPTION "shape4D is an extension for estimating continuous evolution from a discrete collection of shapes, designed to produce realistic anatomical trajectories.")
  set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/shape4D.png")
  set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/shape4D/Screenshots/1.png")
  # Required in Superbuild.cmake
  set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any
  set(SUPERBUILD_TOPLEVEL_PROJECT inner)
  set(EXTENSION_BUILD_SUBDIRECTORY inner-build) # Required for packaging purposes
endif()

#-----------------------------------------------------------------------------
# Superbuild
#-----------------------------------------------------------------------------

# Set to OFF by default, to avoid changing behavior (backward compatibility).
option(${PROJECT_NAME}_SUPERBUILD "Build ${PROJECT_NAME} and the projects it depends on." OFF)
mark_as_advanced(${PROJECT_NAME}_SUPERBUILD)
if(${PROJECT_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
# Inner-build
#-----------------------------------------------------------------------------
include("${CMAKE_CURRENT_SOURCE_DIR}/Inner.cmake")
