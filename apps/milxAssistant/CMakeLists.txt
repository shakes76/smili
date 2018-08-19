############################
# MILXASSISTANT CMAKE FILE #
############################
cmake_minimum_required(VERSION 3.2.0)

PROJECT(milxAssistant)
cmake_policy(SET CMP0053 NEW)

### Include the directory itself as a path to include directories ###
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#######################################################################################
# This sets the windows build that will need the special winmain@16 call. Qt provides #
# this for us in the qtmain.lib file. Using this cmake code will ensure we have it    #
# linked into our build. Not needed on Unix/OS X/Linux which is why we have the       #
# IF(WIN32) conditional.                                                              #
#######################################################################################
IF (WIN32)
    SET (QT_USE_QTMAIN TRUE)
ENDIF (WIN32)

ADD_DEFINITIONS(-DMILX_DLL) #DLL Export MILX

#########################
# Add Headers into List #
#########################
file(GLOB MILXQT_SRC_RCCS ${SMILI_SOURCE_DIR}/*.qrc)

QT5_ADD_RESOURCES(MILXQT_RCCS ${MILXQT_SRC_RCCS})

#################################################
# Iconing and APP variable for fixup            #
# Needs to be before ADD_EXECUTABLE for Windows #
#################################################
ADD_EXECUTABLE(milxAssistant WIN32 MACOSX_BUNDLE milxAssistant.cpp ${MILXQT_RCCS}) #no MACOSX_BUNDLE since custom one created
IF(VTK_QT_VERSION VERSION_GREATER "4")
  qt5_use_modules(milxAssistant Core Gui Widgets WebEngineWidgets)
ENDIF(VTK_QT_VERSION VERSION_GREATER "4")
TARGET_LINK_LIBRARIES(milxAssistant ${QT_LIBRARIES})

# Installs the executable
INSTALL(TARGETS milxAssistant 
	BUNDLE DESTINATION . COMPONENT Runtime
        RUNTIME DESTINATION bin COMPONENT Runtime)
