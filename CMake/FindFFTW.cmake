# - Find the FFTW library
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   FFTW_FOUND               ... true if fftw is found on the system
#   FFTW_LIBRARIES           ... full path to fftw library
#   FFTW_INCLUDE_DIR         ... fftw include directory

if(NOT FFTW_INCLUDE_DIR OR NOT FFTW_LIBRARIES)
    
    #find libs
   if(WIN32)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")
        set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
    endif()
    find_library(
        FFTW_LIB
        NAMES fftw3 fftw3-3
        HINTS 
        ${FFTW_ROOT}
        ENV FFTW_ROOT
		PATH_SUFFIXES lib
    )
    
    #find includes
    find_path(
        FFTW_INCLUDE_DIR
        NAMES fftw3.h
        HINTS 
        ${FFTW_ROOT}
        ENV FFTW_ROOT
        PATH_SUFFIXES include
    )
    
endif()

set(FFTW_LIBRARIES ${FFTW_LIB} ${FFTWF_LIB})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
FFTW_INCLUDE_DIR FFTW_LIBRARIES)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARIES FFTW_LIB FFTWF_LIB)

