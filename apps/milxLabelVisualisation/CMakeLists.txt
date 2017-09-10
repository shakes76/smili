# CMAKE FILE FOR MILXLABELVISUALISATION
cmake_minimum_required(VERSION 3.9)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(milxLabelVisualisation)
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
    find_package(Qt5Widgets REQUIRED QUIET)
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
set(MILXQT_BINDIR ${SMILI_BINARY_DIR})
ADD_DEFINITIONS(-DMILX_DLL) #DLL Export MILX
ADD_DEFINITIONS(-DMILX_PLUGIN_DLL) #DLL Export MILX
IF(WIN32)
	SET(MILXQT_INSTALL_BINARY_DIR "C:/Program Files/MILX/lib")
ELSE(WIN32)
	SET(MILXQT_INSTALL_BINARY_DIR /usr/local/lib)
ENDIF(WIN32)

#######################################
# Include and link external libraries #
#######################################
include_directories (${MILXQT_INCLUDE_DIRS})
link_directories(${MILXQT_INSTALL_BINARY_DIR})

#########################################
# Add executable and link the libraries #
#########################################
ADD_EXECUTABLE(milxLabelVisualisation MACOSX_BUNDLE milxLabelVisualisation.cpp)
TARGET_LINK_LIBRARIES(milxLabelVisualisation milx-Qt milx-SMILI ${QT_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES})

#######################
# Handle UNIX support #
#######################
IF(UNIX)
  configure_file(milxLabelVisualisation.in "milxLabelVisualisation" @ONLY)
  
  install(FILES "${PROJECT_BINARY_DIR}/milxLabelVisualisation" DESTINATION "../../bin" 
    PERMISSIONS WORLD_READ WORLD_EXECUTE)
ENDIF(UNIX)

##############################################
# Expand the command line wild-cards in MSVC #
##############################################
IF(MSVC)
  SET_TARGET_PROPERTIES(milxLabelVisualisation PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

##########################
# Install the executable #
##########################
INSTALL(TARGETS milxLabelVisualisation 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime
	)
