# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

cmake_minimum_required(VERSION 2.6)

PROJECT(finitetransform CXX)

FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtWebkit)
INCLUDE(${QT_USE_FILE})

set(SMILI_DIR ${SMILI_SOURCE_DIR})

# Add Headers into List
file(GLOB SMILI_HEADERS ${SMILI_DIR}/include/*.h)

include_directories (${MILXQT_INCLUDE_DIRS})
link_directories(${SMILI_BINARY_DIR}/lib)

ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
#ADD_DEFINITIONS(-DQT_NO_DEBUG)
ADD_DEFINITIONS(-DQT_SHARED)
ADD_DEFINITIONS(-DMILXQT_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_MAKEDLL) #DLL Plugin Export SMILI
#Win32 requirement for VTK classes in DLLs
#~ ADD_DEFINITIONS (-DvtkFiltering_EXPORTS)

#NTTW Options
#Add compiler flags
# 64-bit Build?
IF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment
  set (NTTW_64 ON) #True
  message("Detected 64-bit system. Using appropriate 64-bit integers.")
ELSEIF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment
  set (NTTW_64 OFF) #False
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8) #64-bit environment

IF(WIN32)
    IF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) #Win64 environment
        set(NTTW_DIR $ENV{ProgramW6432}/NTTW CACHE STRING "NTTW Base Directory")
    ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) # else use Win32 environment
        set(NTTW_DIR $ENV{PROGRAMFILES}/NTTW CACHE STRING "NTTW Base Directory")
    ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432})
ENDIF(WIN32)
IF(UNIX)
    set(NTTW_DIR /usr/local CACHE PATH "NTTW Base Directory")
ENDIF(UNIX)

add_definitions (-DNTTW_DLL)
if (NTTW_64)
	add_definitions (-DNTTW_64)
endif (NTTW_64)
MARK_AS_ADVANCED(NTTW_DIR)
message("Using NTTW from ${NTTW_DIR}")
# Add Headers into List
include_directories (${NTTW_DIR}/include)
link_directories(${NTTW_DIR}/lib)

# this moc's the above variable and appends to the cxx sources
QT4_WRAP_CPP(SSM_MOCS milxQtFiniteTransformPlugin.h)

ADD_LIBRARY(finitetransform SHARED milxQtFiniteTransformPlugin.cpp ${SSM_MOCS})
TARGET_LINK_LIBRARIES(finitetransform nttw milx-SMILI milx-Qt ${QT_LIBRARIES})

install(TARGETS finitetransform
        DESTINATION "plugins"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
