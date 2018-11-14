# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

cmake_minimum_required(VERSION 2.6)

PROJECT(DICOMPlugin CXX)

include(find)

if("${ITK_VERSION_MAJOR}" LESS 4)
    message("DICOM Plugin requires ITK 4 or above for DICOM RT support. Disabling structure set support.")
endif("${ITK_VERSION_MAJOR}" LESS 4)

#~ FIND_PACKAGE(VTK REQUIRED)
#~ INCLUDE(${VTK_USE_FILE})

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

# this moc's the above variable and appends to the cxx sources
QT4_WRAP_CPP(SSM_MOCS milxQtDICOMPlugin.h)

ADD_LIBRARY(DICOMPlugin SHARED milxQtDICOMPlugin.cpp ${SSM_MOCS})
TARGET_LINK_LIBRARIES(DICOMPlugin milx-SMILI milx-Qt ${QT_LIBRARIES})

if(NOT APPLE)
	install(TARGETS DICOMPlugin
	        DESTINATION "plugins"
	        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif()