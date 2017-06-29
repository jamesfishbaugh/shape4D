
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

    if(WIN32)
        # Assume that Slicer is always built in 64bits on Windows
        set(ftp://ftp.fftw.org/pub/fftw/fftw-3.3.5-dll64.zip)
    else()
        set(DOWNLOAD_URL http://fftw.org/fftw-3.3.6-pl2.tar.gz)
    endif()

    set(FFTW_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/FFTW)
    set(FFTW_build ${CMAKE_CURRENT_BINARY_DIR}/FFTW-build)
    include(ExternalProjectForNonCMakeProject)

    set(tcl_SOURCE_DIR FFTW)
    set(${proj}_BUILD_IN_SOURCE 1)

    # environment
    set(_env_script ${CMAKE_BINARY_DIR}/${proj}_Env.cmake)
    ExternalProject_Write_SetBuildEnv_Commands(${_env_script})

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
    set(FFTW_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${_configure_script})

    # build step
    set(_build_script ${CMAKE_BINARY_DIR}/${proj}_build_step.cmake)
    file(WRITE ${_build_script}
"include(\"${_env_script}\")
set(${proj}_WORKING_DIR \"${FFTW_SOURCE_DIR}\")
ExternalProject_Execute(${proj} \"build\" make)
")
    set(FFTW_BUILD_COMMAND ${CMAKE_COMMAND} -P ${_build_script})

    # install step
    set(_install_script ${CMAKE_BINARY_DIR}/${proj}_install_step.cmake)
    file(WRITE ${_install_script}
"include(\"${_env_script}\")
set(${proj}_WORKING_DIR \"${FFTW_SOURCE_DIR}\")
ExternalProject_Execute(${proj} \"install\" make install)
")
    set(FFTW_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${_install_script})


    ExternalProject_Add(${proj}
        ${${proj}_EP_ARGS}
        URL ${DOWNLOAD_URL}
        UPDATE_COMMAND "" # Disable update
        SOURCE_DIR ${FFTW_SOURCE_DIR}
        BUILD_IN_SOURCE ${FFTW_BUILD_IN_SOURCE}
        CONFIGURE_COMMAND ${FFTW_CONFIGURE_COMMAND}
        BUILD_COMMAND ${FFTW_BUILD_COMMAND}
        INSTALL_COMMAND ${FFTW_INSTALL_COMMAND}
        DEPENDS
          ${${proj}_DEPENDS}
    )
    set(${proj}_ROOT ${CMAKE_BINARY_DIR}/${proj}-build)
    set(INSTALL_LIBRARIES
        ${${proj}_ROOT}/lib/libfftw3.so
        ${${proj}_ROOT}/lib/libfftw3.so.3
        ${${proj}_ROOT}/lib/libfftw3.so.3.5.6
    )
else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDS})
endif()

mark_as_superbuild(${proj}_ROOT:PATH)
