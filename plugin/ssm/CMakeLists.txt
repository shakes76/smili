# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

cmake_minimum_required(VERSION 2.6)

PROJECT(SSMPlugin CXX)

include(find)

FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtWebkit)
INCLUDE(${QT_USE_FILE})

set(SMILI_DIR ${SMILI_SOURCE_DIR})

# Add Headers into List
file(GLOB SMILI_HEADERS ${SMILI_DIR}/include/*.h)

include_directories (${MILXQT_INCLUDE_DIRS})
link_directories(${SMILI_BINARY_DIR}/lib)

#MILX-MSK Options
IF(WIN32)
    IF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) #Win64 environment
        set(MILX-MSK_DIR $ENV{ProgramW6432}/MILX_MSK CACHE STRING "MILX-MSK Base Directory")
    ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) # else use Win32 environment
        set(MILX-MSK_DIR $ENV{PROGRAMFILES}/MILX_MSK CACHE STRING "MILX-MSK Base Directory")
    ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432})
ENDIF(WIN32)
IF(UNIX)
    set(MILX-MSK_DIR /usr/local CACHE PATH "MILX-MSK Base Directory")
ENDIF(UNIX)
MARK_AS_ADVANCED(MILX-MSK_DIR)
message("Using MSK from ${MILX-MSK_DIR}")
# Add Headers into List
file(GLOB MILX-MSK_HEADERS ${MILX-MSK_DIR}/include/*.h)
include_directories (${MILX-MSK_DIR}/include)
include_directories (${MILX-MSK_DIR}/include/milxSSM3D) #for SSMs and RSSMs
include_directories (${MILX-MSK_DIR}/include/itk-ext/Modules/IO) # for INR format
include_directories (${MILX-MSK_DIR}/include/vtk-ext/Modules/Hybrid) # for WPCA
link_directories(${MILX-MSK_DIR}/lib)

set(MILX-VIEW_DIR /usr/local/milxview CACHE PATH "MilxView  Base Directory")
message("Using MilxView from ${MILX-VIEW_DIR}")
MARK_AS_ADVANCED(MILX-VIEW_DIR)
link_directories(${MILX-VIEW_DIR}/lib) #For milx-itkIO, milx-vtkHybrid

ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
#ADD_DEFINITIONS(-DQT_NO_DEBUG)
ADD_DEFINITIONS(-DQT_SHARED)
ADD_DEFINITIONS(-DMILXQT_DLL) #DLL Export milxQt
ADD_DEFINITIONS(-DMILXQT_PLUGIN_DLL) #DLL Export milxQt
ADD_DEFINITIONS(-DMILXQT_PLUGIN_MAKEDLL) #DLL Plugin Export milxQt

# this moc's the above variable and appends to the cxx sources
QT4_WRAP_CPP(SSM_MOCS milxQtSSMPlugin.h milxQtShapeModel.h milxQtRobustShapeModel.h)


set(SSM_USE_LOCAL_FILE false CACHE BOOL "Use a local version of the vtkPCAAnalysisTimepointFilter*, vtkWeightedPCAAnalysisFilter* etc. files")
IF (SSM_USE_LOCAL_FILE)
    set(SSM_LOCAL_FILES vtkPCAAnalysisTimepointFilter.cxx vtkWeightedPCAAnalysisFilter.cxx)
ELSE (SSM_LOCAL_FILES)
    set(SSM_LOCAL_FILES "")
ENDIF (SSM_USE_LOCAL_FILE)

ADD_LIBRARY(SSMPlugin SHARED milxQtShapeModel.cpp milxQtRobustShapeModel.cpp milxQtSSMPlugin.cpp ${SSM_LOCAL_FILES} ${SSM_MOCS})

IF (SSM_USE_LOCAL_FILE)
  TARGET_LINK_LIBRARIES(SSMPlugin vtk-ext milx-Qt ${QT_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${ITK_REVIEW_LIBRARIES})
ELSE (SSM_USE_LOCAL_FILE)
  TARGET_LINK_LIBRARIES(SSMPlugin milx-vtkHybrid vtk-ext milx-Qt ${QT_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${ITK_REVIEW_LIBRARIES})
ENDIF (SSM_USE_LOCAL_FILE)

if(BUILD_SSM_PLUGIN)
  install(TARGETS SSMPlugin
          DESTINATION "plugins"
          PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif(BUILD_SSM_PLUGIN)
