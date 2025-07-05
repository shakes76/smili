cmake_minimum_required(VERSION 3.12)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(PythonPlugin CXX)

#~ set(Python3_ROOT_DIR "D:/WPy64-310111-JAX/python-3.10.11.amd64")
set(Python3_ROOT_DIR "$ENV{PYTHON_PATH}")
MARK_AS_ADVANCED(Python3_ROOT_DIR)

find_package(Python3 COMPONENTS Development)
MESSAGE("Found Python 3 at", ${Python3_ROOT_DIR})

set(MILXQT_DIR ${SMILI_SOURCE_DIR})

#PythonQt Options
IF(WIN32)
    IF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) #Win64 environment
        set(PYTHONQT_DIR $ENV{ProgramW6432}/PythonQt CACHE PATH "Python Qt Base Directory")
        set(PYTHONQT_INCLUDE_PATH $ENV{ProgramW6432}/PythonQt/include CACHE PATH "Python Qt Include Directory")
    ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) # else use Win32 environment
        set(PYTHONQT_DIR $ENV{PROGRAMFILES}/PythonQt CACHE PATH "Python Qt Base Directory")
        set(PYTHONQT_INCLUDE_PATH $ENV{PROGRAMFILES}/PythonQt/include CACHE PATH "Python Qt Include Directory")
    ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432})
ENDIF(WIN32)
IF(UNIX)
	set(PYTHONQT_DIR /usr CACHE STRING "Python Qt Base Directory")
	set(PYTHONQT_INCLUDE_PATH /usr/include/PythonQt CACHE STRING "Python Qt Base Directory") #Ubuntu 12.04
ENDIF(UNIX)
MARK_AS_ADVANCED(PYTHONQT_DIR)
MARK_AS_ADVANCED(PYTHONQT_INCLUDE_PATH)
include_directories (${PYTHONQT_INCLUDE_PATH})
link_directories(${PYTHONQT_DIR}/lib)
file(GLOB PYTHONQT_LIBS ${PYTHONQT_DIR}/lib/*${CMAKE_STATIC_LIBRARY_SUFFIX}*)
message("Using Python Qt from ${PYTHONQT_DIR}")

# Add Headers into List
file(GLOB MILXQT_HEADERS ${MILXQT_DIR}/include/*.h)

include_directories (${MILXQT_INCLUDE_DIRS} ${Python3_INCLUDE_DIRS} ${PYTHONQT_INCLUDE_PATH})
link_directories(${MILXQT_BINARY_DIR}/lib ${PYTHONQT_DIR}/lib ${Python3_LIBRARY_DIRS})

ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
#ADD_DEFINITIONS(-DQT_NO_DEBUG)
ADD_DEFINITIONS(-DQT_SHARED)
ADD_DEFINITIONS(-DMILXQT_DLL) #DLL Export MILX
ADD_DEFINITIONS(-DMILXQT_PLUGIN_DLL) #DLL Export MILXQT
ADD_DEFINITIONS(-DMILXQT_PLUGIN_MAKEDLL) #DLL Plugin Export MILXQT

# this moc's the above variable and appends to the cxx sources
IF(VTK_QT_VERSION VERSION_GREATER "4")
    QT5_WRAP_CPP(SSM_MOCS milxQtPythonPlugin.h milxQtPythonConsole.h)
ELSEIF(VTK_QT_VERSION VERSION_GREATER "4")
    QT4_WRAP_CPP(SSM_MOCS milxQtPythonPlugin.h milxQtPythonConsole.h)
ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

ADD_LIBRARY(PythonPlugin SHARED milxQtPythonConsole.cpp milxQtPythonPlugin.cpp ${SSM_MOCS} ${PYTHON_LIBRARIES})
TARGET_LINK_LIBRARIES(PythonPlugin milx-Qt ${QT_LIBRARIES} ${Python3_LIBRARIES} ${PYTHONQT_LIBS})

if(BUILD_PYTHON_PLUGIN)
  install(TARGETS PythonPlugin
          DESTINATION "plugins"
          PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif(BUILD_PYTHON_PLUGIN)
