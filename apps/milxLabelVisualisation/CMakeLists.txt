# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(milxLabelVisualisation)

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
