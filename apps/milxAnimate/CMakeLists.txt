# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(milxAnimate)

set(MILXQT_DIR ${SMILI_SOURCE_DIR})
set(MILXQT_BINDIR ${SMILI_BINARY_DIR})
ADD_DEFINITIONS(-DMILX_DLL) #DLL Export MILX
ADD_DEFINITIONS(-DMILX_PLUGIN_DLL) #DLL Export MILX
IF(WIN32)
	SET(MILXQT_INSTALL_BINARY_DIR "C:/Program Files/MILX/lib")
ELSE(WIN32)
	SET(MILXQT_INSTALL_BINARY_DIR /usr/local/lib)
ENDIF(WIN32)

#MILX-MSK Options
IF(WIN32)
	set(MILX-MSK_DIR $ENV{PROGRAMFILES}/MILX_MSK)
#~     message("${MILX-MSK_DIR}")
ENDIF(WIN32)
IF(UNIX)
	set(MILX-MSK_DIR /usr/local)
ENDIF(UNIX)

# Add Headers into List
file(GLOB MILXQT_HEADERS ${MILXQT_DIR}/include/*.h)
file(GLOB MILXQT_LIBS ${SMILI_BINARY_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB MILXQT_PLUGINS ${SMILI_PLUGINS_PATH}/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB PYTHONQT_LIBS ${PYTHONQT_DIR}/lib/*PythonQt*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB VTK_QT_LIBS ${VTK_LIBRARY_DIRS}/*Q*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB MSK_LIBS ${MILX-MSK_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)

include_directories (${MILXQT_INCLUDE_DIRS} ${MILXQT_DIR}/plugin/animate)
link_directories(${MILXQT_INSTALL_BINARY_DIR} ${SMILI_PLUGINS_PATH})

ADD_EXECUTABLE(milxAnimate MACOSX_BUNDLE milxAnimate.cpp)
TARGET_LINK_LIBRARIES(milxAnimate milx-Qt milx-SMILI ${QT_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES} animatePlugin)

# expand the command line wild-cards in msvc
IF(MSVC)
  SET_TARGET_PROPERTIES(milxAnimate PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

MESSAGE("Animation: milxAnimate app is being built. Please ensure the animate plugin is built and linker can locate it.")
MESSAGE("Animation: For linux, setup LD_LIBRARY_PATH or ldconfig appropriately.")

# Installs the executable
INSTALL(TARGETS milxAnimate 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime
	)
