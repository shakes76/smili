#
# This file is part of SMILI.
#
set(PYTHON_BINDING_NAMESPACE "PySMILI")

set(Python3_ROOT_DIR "D:/WPy64-31090/python-3.10.9.amd64") # Replace with the actual path, such as your venv
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
# Use provided python interpreter if given.
if(NOT Python3_EXECUTABLE)
    if(WIN32 AND "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        find_program(Python3_EXECUTABLE "python_d")
        if(NOT Python3_EXECUTABLE)
            message(FATAL_ERROR
                "A debug Python interpreter could not be found, which is a requirement when "
                "building this example in a debug configuration. Make sure python_d.exe is in "
                "PATH.")
        endif()
    else()
        find_program(Python3_EXECUTABLE "python")
        if(NOT Python3_EXECUTABLEr)
            message(FATAL_ERROR
                "No Python interpreter could be found. Make sure python is in PATH.")
        endif()
    endif()
endif()
message("Using Python ${Python3_EXECUTABLE}")
cmake_path(GET Python3_EXECUTABLE PARENT_PATH Python3_BINPATH)
message("Using Python path ${Python3_BINPATH}")
message("Using Python stdlib ${Python3_STDLIB}")
message("Using Python sitelib ${Python3_SITELIB}")

# Just to fix warnings with --warn-uninitialized
if(NOT DEFINED SHIBOKEN_CUSTOM_PREFIX) #look for shiboken in a custom location
    set(SHIBOKEN_CUSTOM_PREFIX "")
endif()
if(NOT DEFINED PYSIDE_CUSTOM_PREFIX) #look for pyside in a custom location
    set(PYSIDE_CUSTOM_PREFIX "")
endif()

if(${PROJECT_NAME}_QT6)
    set(PYSIDE_MAJOR_VERSION "6")
    #set(PYTHON_BINDING_NAMESPACE "${PYTHON_BINDING_NAMESPACE}Qt${PYSIDE_MAJOR_VERSION}")
    set(QtWidgets_VERSION ${Qt6Widgets_VERSION})
else()
    set(PYSIDE_MAJOR_VERSION "2")
    set(QtWidgets_VERSION ${Qt5Widgets_VERSION})
endif()
message("Using PySide ${PYSIDE_MAJOR_VERSION}")

find_package(PySide${PYSIDE_MAJOR_VERSION} ${QtWidgets_VERSION} EXACT REQUIRED)
find_package(Shiboken${PYSIDE_MAJOR_VERSION} REQUIRED)

if(NOT WIN32)
    set(Python3_VERSION_MAJORMINOR "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")
    set(BINDINGS_DIR "python${Python3_VERSION_MAJORMINOR}/${PYTHON_BINDING_NAMESPACE}")
else()
    set(BINDINGS_DIR "${PYTHON_BINDING_NAMESPACE}")
endif()
# set(${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX "${${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX}/${BINDINGS_DIR}")
#Convert a native path stored to a CMake-style path
file(TO_CMAKE_PATH "${Python3_SITELIB}/${BINDINGS_DIR}" ${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX)
#set(${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX "${Python3_SITELIB}/${BINDINGS_DIR}")

#if(NOT ${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX)
#    # cmake-lint: disable=C0103
#    set(${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX
#        ${INSTALL_LIBRARY_DIR}
#        CACHE FILEPATH "Custom path to install python bindings."
#    )
#endif()

# Enable rpaths so that the built shared libraries find their dependencies.
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH ${SHIBOKEN_BASEDIR} ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(PySide${PYSIDE_MAJOR_VERSION}ModuleBuild)
add_subdirectory(PySMILI)
#~ if(${PROJECT_NAME}_TESTS)
#~   add_subdirectory(tests)
#~ endif()
