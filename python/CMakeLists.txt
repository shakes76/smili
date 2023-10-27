#
# This file is part of SMILI.
#
set(PYTHON_BINDING_NAMESPACE "PySMILI")

#~ find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
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

# Just to fix warnings with --warn-uninitialized
if(NOT DEFINED SHIBOKEN_CUSTOM_PREFIX) #look for shiboken in a custom location
  set(SHIBOKEN_CUSTOM_PREFIX "")
endif()
if(NOT DEFINED PYSIDE_CUSTOM_PREFIX) #look for pyside in a custom location
  set(PYSIDE_CUSTOM_PREFIX "")
endif()

if(${PROJECT_NAME}_QT6)
  set(PYSIDE_MAJOR_VERSION "6")
  set(PYTHON_BINDING_NAMESPACE "${PYTHON_BINDING_NAMESPACE}Qt${PYSIDE_MAJOR_VERSION}")
  set(QtWidgets_VERSION ${Qt6Widgets_VERSION})
  message("Using PySide ${PYSIDE_MAJOR_VERSION}")
else()
  set(PYSIDE_MAJOR_VERSION "2")
  set(QtWidgets_VERSION ${Qt5Widgets_VERSION})
endif()

#~ set(Python_ROOT_DIR "C:/Users/uqscha22/.conda/envs/vtk")
find_package(PySide${PYSIDE_MAJOR_VERSION} ${QtWidgets_VERSION} EXACT REQUIRED)
find_package(Shiboken${PYSIDE_MAJOR_VERSION} REQUIRED)

#~ if(NOT PYTHON_BINDINGS_INSTALL_PREFIX)
#~   set(PYTHON_BINDINGS_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
#~ endif()
#~ set(BINDINGS_DIR "${INSTALL_LIBRARY_DIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/${PYTHON_BINDING_NAMESPACE}")
set(BINDINGS_DIR "${Python3_BINPATH}/../lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/${PYTHON_BINDING_NAMESPACE}")
set(PYTHON_BINDINGS_INSTALL_PREFIX "${BINDINGS_DIR}")

include(PySide${PYSIDE_MAJOR_VERSION}ModuleBuild)
add_subdirectory(PySMILI)
#~ if(${PROJECT_NAME}_TESTS)
#~   add_subdirectory(tests)
#~ endif()
