# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

cmake_minimum_required(VERSION 2.6)

PROJECT(PythonPlugin CXX)

FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtWebkit)
INCLUDE(${QT_USE_FILE})

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
message("Using Python Qt from ${PYTHONQT_DIR}")

# Add Headers into List
file(GLOB MILXQT_HEADERS ${MILXQT_DIR}/include/*.h)

include_directories (${MILXQT_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${PYTHONQT_INCLUDE_PATH})
link_directories(${MILXQT_BINARY_DIR}/lib ${PYTHONQT_DIR}/lib)

ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
#ADD_DEFINITIONS(-DQT_NO_DEBUG)
ADD_DEFINITIONS(-DQT_SHARED)
ADD_DEFINITIONS(-DMILXQT_DLL) #DLL Export MILX
ADD_DEFINITIONS(-DMILXQT_PLUGIN_DLL) #DLL Export MILXQT
ADD_DEFINITIONS(-DMILXQT_PLUGIN_MAKEDLL) #DLL Plugin Export MILXQT

# this moc's the above variable and appends to the cxx sources
QT4_WRAP_CPP(SSM_MOCS milxQtPythonPlugin.h milxQtPythonConsole.h)

ADD_LIBRARY(PythonPlugin SHARED milxQtPythonConsole.cpp milxQtPythonPlugin.cpp ${SSM_MOCS} ${PYTHON_LIBRARIES})
TARGET_LINK_LIBRARIES(PythonPlugin milx-Qt ${QT_LIBRARIES} PythonQt ${PYTHON_LIBRARIES})

if(BUILD_PYTHON_PLUGIN)
  install(TARGETS PythonPlugin
          DESTINATION "plugins"
          PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif(BUILD_PYTHON_PLUGIN)
