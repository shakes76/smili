# CMAKE FILE FOR MILX ASSISTANT
cmake_minimum_required(VERSION 3.9)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(milxAssistant)

#######################################################################################
# This sets the windows build that will need the special winmain@16 call. Qt provides #
#  this for us in the qtmain.lib file. Using this cmake code will ensure we have it   #
#  linked into our build. Not needed on Unix/OS X/Linux which is why we have the      #
#  IF(WIN32) conditional.                                                             #
#######################################################################################
IF (WIN32)
    SET (QT_USE_QTMAIN TRUE)
ENDIF (WIN32)

###############################
# Find all necessary packages #
###############################
FIND_PACKAGE(Qt5 REQUIRED Widgets WebEngineWidgets)
ADD_DEFINITIONS(-DMILX_DLL) #DLL Export MILX

#########################
# Add Headers into List #
#########################
file(GLOB MILXQT_SRC_RCCS ${SMILI_SOURCE_DIR}/*.qrc)
QT5_ADD_RESOURCES(MILXQT_RCCS ${MILXQT_SRC_RCCS})

#####################################################################################
# Iconing and APP variable for fixup, needs to be before ADD_EXECUTABLE for Windows #
#####################################################################################
ADD_EXECUTABLE(milxAssistant WIN32 MACOSX_BUNDLE milxAssistant.cpp ${MILXQT_RCCS}) # No MACOSX_BUNDLE since custom one created
TARGET_LINK_LIBRARIES(milxAssistant Qt5::Widgets Qt5::WebEngineWidgets)

##########################
# Install the executable #
##########################
INSTALL(TARGETS milxAssistant 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime)
