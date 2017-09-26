
set(proj FFTW)

# Set dependency list
set(${proj}_DEPENDS "")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  message(FATAL_ERROR "Enabling ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj} is not supported !")
endif()

# Sanity checks
if(DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR})
  message(FATAL_ERROR "FFTW_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED ${proj}_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

    if(NOT DEFINED git_protocol)
        set(git_protocol "git")
    endif()

    set(${proj}_CONFIGURE_COMMAND "")
    set(${proj}_INSTALL_COMMAND "")

    set(${proj}_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
    set(${proj}_ROOT ${${proj}_SOURCE_DIR})
    
    # Needed to generate and run the configure/build/install scripts
    include(ExternalProjectForNonCMakeProject RESULT_VARIABLE ExternalProjectForNonCMakeProject_path)

    # environment
    set(_env_script ${CMAKE_BINARY_DIR}/${proj}_Env.cmake)
    ExternalProject_Write_SetBuildEnv_Commands(${_env_script})

    if(WIN32) #Windows
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
		set(${proj}_WORKING_DIR \"${${proj}_SOURCE_DIR}\")
		ExternalProject_Execute(${proj} \"build\" \"${_vs_bin_path}/lib.exe\" /machine:x64 /def:libfftw3-3.def)
	")
        set(${proj}_BUILD_COMMAND ${CMAKE_COMMAND} -P ${_build_script})
        
        set(${proj}_INSTALL_LIBRARIES
            ${${proj}_ROOT}/libfftw3-3.dll
            ${${proj}_ROOT}/libfftw3-3.lib
        )

    else() # Linux or MAC
        set(${proj}_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj})
        set(${proj}_build ${CMAKE_CURRENT_BINARY_DIR}/${proj}-build)
        set(DOWNLOAD_URL http://fftw.org/fftw-3.3.6-pl2.tar.gz)

        # configure step
        set(_configure_script ${CMAKE_BINARY_DIR}/${proj}_configure_step.cmake)
        file(WRITE ${_configure_script}
		"include(${ExternalProjectForNonCMakeProject_path})
		set(CMAKE_BINARY_DIR ${CMAKE_BINARY_DIR})
		set(${proj}_WORKING_DIR \"${${proj}_SOURCE_DIR}\")
		ExternalProject_Execute(${proj} \"configure\" sh configure
		  --prefix=${${proj}_build}
		  --enable-shared --enable-static=no
		  )
	")

        # build step
        set(_build_script ${CMAKE_BINARY_DIR}/${proj}_build_step.cmake)
        file(WRITE ${_build_script}
		"include(\"${_env_script}\")
		set(${proj}_WORKING_DIR \"${${proj}_SOURCE_DIR}\")
		ExternalProject_Execute(${proj} \"build\" make)
	")

        # install step
        set(_install_script ${CMAKE_BINARY_DIR}/${proj}_install_step.cmake)
        file(WRITE ${_install_script}
		"include(\"${_env_script}\")
		set(${proj}_WORKING_DIR \"${FFTW_SOURCE_DIR}\")
		ExternalProject_Execute(${proj} \"install\" make install)
	")

        set(${proj}_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${_configure_script})
        set(${proj}_BUILD_COMMAND ${CMAKE_COMMAND} -P ${_build_script})
        set(${proj}_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${_install_script})


	if(APPLE)

        	set(${proj}_INSTALL_LIBRARIES
            	  ${${proj}_ROOT}/lib/libfftw3.3.dylib
            	  ${${proj}_ROOT}/lib/libfftw3.dylib
            	  ${${proj}_ROOT}/lib/libfftw3.la
        	)

		set(${proj}_LIB
	    	  ${${proj}_ROOT}/lib/libfftw3.dylib)

	else()
        	set(${proj}_INSTALL_LIBRARIES
            	  ${${proj}_ROOT}/lib/libfftw3.so
            	  ${${proj}_ROOT}/lib/libfftw3.so.3
            	  ${${proj}_ROOT}/lib/libfftw3.so.3.5.6
        	)

		set(${proj}_LIB
	    	  ${${proj}_ROOT}/lib/libfftw3.so)

	endif()

	set(${proj}_INCLUDE_DIR 
	    ${${proj}_ROOT}/include)
    endif()
    
    set(${proj}_BUILD_IN_SOURCE 1)

    ExternalProject_Add(${proj}
        ${${proj}_EP_ARGS}
        URL ${DOWNLOAD_URL}
        UPDATE_COMMAND "" # Disable update
        SOURCE_DIR ${${proj}_SOURCE_DIR}
        BUILD_IN_SOURCE ${${proj}_BUILD_IN_SOURCE}
        CONFIGURE_COMMAND "${${proj}_CONFIGURE_COMMAND}"
        BUILD_COMMAND "${${proj}_BUILD_COMMAND}"
        INSTALL_COMMAND "${${proj}_INSTALL_COMMAND}"
        DEPENDS
          ${${proj}_DEPENDS}
    )
else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDS})
endif()

mark_as_superbuild(${proj}_ROOT:PATH)
mark_as_superbuild(${proj}_INSTALL_LIBRARIES:STRING)
mark_as_superbuild(${proj}_INCLUDE_DIR:STRING)
mark_as_superbuild(${proj}_LIB:STRING)
