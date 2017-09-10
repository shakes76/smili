# CMAKE FILE FOR MILXOVERLAY
cmake_minimum_required(VERSION 3.9)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(milxOverlay)
cmake_policy(SET CMP0020 NEW)

#############################################
# Find VTK and ITK packages and directories #
#############################################
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

#######################################################################################
# This sets the windows build that will need the special winmain@16 call. Qt provides #
# this for us in the qtmain.lib file. Using this cmake code will ensure we have it    #
# linked into our build. Not needed on Unix/OS X/Linux which is why we have the       #
# IF(WIN32) conditional. 							      #
#######################################################################################
IF (WIN32)
    SET (QT_USE_QTMAIN TRUE)
ENDIF (WIN32)

####################################################################
# Detect Qt version (depends on VTK-Qt) and get required libraries #
####################################################################
IF(VTK_QT_VERSION VERSION_GREATER "4")
    FIND_PACKAGE(Qt5Widgets REQUIRED QUIET)
    include_directories(${Qt5Widgets_INCLUDE_DIRS})
    add_definitions(${Qt5Widgets_DEFINITIONS})
    qt5_wrap_cpp(QVTKLibMocSrcs ${QVTKMocHeaders})
    set(QT_LIBRARIES ${Qt5Widgets_LIBRARIES})
ELSE(VTK_QT_VERSION VERSION_GREATER "4")
    message("VTK Qt Version: ${VTK_QT_VERSION}")
ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

###############
# Set up MILX #
###############
set(MILXQT_DIR ${SMILI_SOURCE_DIR})
ADD_DEFINITIONS(-DMILX_DLL) #DLL Export MILX
ADD_DEFINITIONS(-DMILX_PLUGIN_DLL) #DLL Export MILX
IF(WIN32)
	SET(MILXQT_INSTALL_BINARY_DIR "C:/Program Files/MILX/lib")
ELSE(WIN32)
	SET(MILXQT_INSTALL_BINARY_DIR /usr/local/lib)
ENDIF(WIN32)

####################
# MILX-MSK Options #
####################
IF(WIN32)
	set(MILX-MSK_DIR $ENV{PROGRAMFILES}/MILX_MSK)
#~     message("${MILX-MSK_DIR}")
ENDIF(WIN32)
IF(UNIX)
	set(MILX-MSK_DIR /usr/local)
ENDIF(UNIX)

#########################
# Add Headers into List #
#########################
file(GLOB MILXQT_HEADERS ${MILXQT_DIR}/include/*.h)
file(GLOB MILXQT_LIBS ${SMILI_BINARY_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB MILXQT_PLUGINS ${SMILI_BINARY_DIR}/plugins/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB PYTHONQT_LIBS ${PYTHONQT_DIR}/lib/*PythonQt*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB VTK_QT_LIBS ${VTK_LIBRARY_DIRS}/*Q*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB MSK_LIBS ${MILX-MSK_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)

#########################################
# Include and link external directories #
#########################################
include_directories (${MILXQT_INCLUDE_DIRS})
link_directories(${MILXQT_INSTALL_BINARY_DIR})

#########################################
# Add executable and link the libraries #
#########################################
ADD_EXECUTABLE(milxOverlay MACOSX_BUNDLE milxOverlay.cpp)
TARGET_LINK_LIBRARIES(milxOverlay milx-Qt milx-SMILI ${QT_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES})

#######################
# Handle UNIX support #
#######################
IF(UNIX)
  configure_file(milxOverlay.in "milxOverlay" @ONLY)
  
  install(FILES "${PROJECT_BINARY_DIR}/milxOverlay" DESTINATION "../../bin" 
    PERMISSIONS WORLD_READ WORLD_EXECUTE)
ENDIF(UNIX)

##############################################
# Expand the command line wild-cards in MSVC #
##############################################
IF(MSVC)
  SET_TARGET_PROPERTIES(milxOverlay PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

##########################
# Install the executable #
##########################
INSTALL(TARGETS milxOverlay 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime
	)
