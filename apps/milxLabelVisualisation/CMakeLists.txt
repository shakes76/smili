cmake_minimum_required(VERSION 2.8)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(milxLabelVisualisation)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

#--------------------------------------------------------------------------------
# This sets the windows build that will need the special winmain@16 call. Qt provides
#  this for us in the qtmain.lib file. Using this cmake code will ensure we have it
#  linked into our build. Not needed on Unix/OS X/Linux which is why we have the 
#  IF(WIN32) conditional.
IF (WIN32)
    SET (QT_USE_QTMAIN TRUE)
ENDIF (WIN32)

set(MILXQT_DIR ${SMILI_SOURCE_DIR})
set(MILXQT_BINDIR ${SMILI_BINARY_DIR})
ADD_DEFINITIONS(-DMILX_DLL) #DLL Export MILX
ADD_DEFINITIONS(-DMILX_PLUGIN_DLL) #DLL Export MILX
IF(WIN32)
	SET(MILXQT_INSTALL_BINARY_DIR "C:/Program Files/MILX/lib")
ELSE(WIN32)
	SET(MILXQT_INSTALL_BINARY_DIR /usr/local/lib)
ENDIF(WIN32)

include_directories (${MILXQT_INCLUDE_DIRS})
link_directories(${MILXQT_INSTALL_BINARY_DIR})

ADD_EXECUTABLE(milxLabelVisualisation MACOSX_BUNDLE milxLabelVisualisation.cpp)
TARGET_LINK_LIBRARIES(milxLabelVisualisation milx-Qt milx-SMILI ${QT_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES})

IF(UNIX)
  configure_file(milxLabelVisualisation.in "milxLabelVisualisation" @ONLY)
  
  install(FILES "${PROJECT_BINARY_DIR}/milxLabelVisualisation" DESTINATION "../../bin" 
    PERMISSIONS WORLD_READ WORLD_EXECUTE)
ENDIF(UNIX)

# expand the command line wild-cards in msvc
IF(MSVC)
  SET_TARGET_PROPERTIES(milxLabelVisualisation PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

# Installs the executable
INSTALL(TARGETS milxLabelVisualisation 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime
	)
