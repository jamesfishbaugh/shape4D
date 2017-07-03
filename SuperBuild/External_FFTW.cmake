
set(proj FFTW)

# Set dependency list
set(${proj}_DEPENDS "")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  message(FATAL_ERROR "Enabling ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj} is not supported !")
endif()

# Sanity checks
if(DEFINED FFTW_DIR AND NOT EXISTS ${FFTW_DIR})
  message(FATAL_ERROR "FFTW_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED ${proj}_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

    if(NOT DEFINED git_protocol)
        set(git_protocol "git")
    endif()

    set(FFTW_CONFIGURE_COMMAND "")
    set(FFTW_INSTALL_COMMAND "")

    set(FFTW_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
    set(${proj}_ROOT ${FFTW_SOURCE_DIR})
    
    # Needed to generate and run the configure/build/install scripts
    include(ExternalProjectForNonCMakeProject)

    # environment
    set(_env_script ${CMAKE_BINARY_DIR}/${proj}_Env.cmake)
    ExternalProject_Write_SetBuildEnv_Commands(${_env_script})
    
    if(WIN32)
        # Assume that Slicer is always built in 64bits on Windows
        set(DOWNLOAD_URL ftp://ftp.fftw.org/pub/fftw/fftw-3.3.5-dll64.zip)
        # We need to create libfftw3-3.lib that is not included in the zip file
        # See https://stackoverflow.com/questions/38095244/\
        # how-to-get-path-to-provided-visual-studio-lib-exe-executable-in-cmake
        get_filename_component(_vs_bin_path "${CMAKE_LINKER}" DIRECTORY)
        # build step
        set(_build_script ${CMAKE_BINARY_DIR}/${proj}_build_step.cmake)
        file(WRITE ${_build_script}
"include(\"${_env_script}\")
set(${proj}_WORKING_DIR \"${FFTW_SOURCE_DIR}\")
ExternalProject_Execute(${proj} \"build\" \"${_vs_bin_path}/lib.exe\" /machine:x64 /def:libfftw3-3.def)
")
        set(FFTW_BUILD_COMMAND ${CMAKE_COMMAND} -P ${_build_script})
        
        set(FFTW_INSTALL_LIBRARIES
            ${${proj}_ROOT}/libfftw3-3.dll
            ${${proj}_ROOT}/libfftw3-3.lib
        )
    elseif(APPLE)
        # See Utilies/README to generate this package
        set(DOWNLOAD_URL "http://slicer.kitware.com/midas3/download/?items=296662,1")
        
        set(FFTW_INSTALL_LIBRARIES
            ${${proj}_ROOT}/lib/libfftw3.3.dylib
            ${${proj}_ROOT}/lib/libfftw3.dylib
            ${${proj}_ROOT}/lib/libfftw3.la
        )
        set(FFTW_BUILD_COMMAND "")
    else()
        set(FFTW_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/FFTW)
        set(FFTW_build ${CMAKE_CURRENT_BINARY_DIR}/FFTW-build)
        set(DOWNLOAD_URL http://fftw.org/fftw-3.3.6-pl2.tar.gz)
    
        # configure step
        set(_configure_script ${CMAKE_BINARY_DIR}/${proj}_configure_step.cmake)
        file(WRITE ${_configure_script}
"include(\"${_env_script}\")
set(${proj}_WORKING_DIR \"${FFTW_SOURCE_DIR}\")
ExternalProject_Execute(${proj} \"configure\" sh configure
  --prefix=${FFTW_build}
  --enable-shared --enable-static=no
  )
")

        # build step
        set(_build_script ${CMAKE_BINARY_DIR}/${proj}_build_step.cmake)
        file(WRITE ${_build_script}
"include(\"${_env_script}\")
set(${proj}_WORKING_DIR \"${FFTW_SOURCE_DIR}\")
ExternalProject_Execute(${proj} \"build\" make)
")

        # install step
        set(_install_script ${CMAKE_BINARY_DIR}/${proj}_install_step.cmake)
        file(WRITE ${_install_script}
"include(\"${_env_script}\")
set(${proj}_WORKING_DIR \"${FFTW_SOURCE_DIR}\")
ExternalProject_Execute(${proj} \"install\" make install)
")

        set(FFTW_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${_configure_script})
        set(FFTW_BUILD_COMMAND ${CMAKE_COMMAND} -P ${_build_script})
        set(FFTW_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${_install_script})

        set(FFTW_INSTALL_LIBRARIES
            ${${proj}_ROOT}/lib/libfftw3.so
            ${${proj}_ROOT}/lib/libfftw3.so.3
            ${${proj}_ROOT}/lib/libfftw3.so.3.5.6
        )
    endif()
    
    set(${proj}_BUILD_IN_SOURCE 1)

    ExternalProject_Add(${proj}
        ${${proj}_EP_ARGS}
        URL ${DOWNLOAD_URL}
        UPDATE_COMMAND "" # Disable update
        SOURCE_DIR ${FFTW_SOURCE_DIR}
        BUILD_IN_SOURCE ${FFTW_BUILD_IN_SOURCE}
        CONFIGURE_COMMAND "${FFTW_CONFIGURE_COMMAND}"
        BUILD_COMMAND "${FFTW_BUILD_COMMAND}"
        INSTALL_COMMAND "${FFTW_INSTALL_COMMAND}"
        DEPENDS
          ${${proj}_DEPENDS}
    )
else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDS})
endif()

mark_as_superbuild(${proj}_ROOT:PATH)
mark_as_superbuild(FFTW_INSTALL_LIBRARIES:STRING)
