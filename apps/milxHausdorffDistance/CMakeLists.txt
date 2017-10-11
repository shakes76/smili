cmake_minimum_required(VERSION 2.8)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(milxHausdorffDistance)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

ADD_DEFINITIONS(-DMILX_DLL) #DLL Export MILX
ADD_DEFINITIONS(-DMILX_PLUGIN_DLL) #DLL Export MILX

include_directories (${SMILI_INCLUDE_DIRS})
link_directories(${MILXQT_INSTALL_BINARY_DIR})

ADD_EXECUTABLE(milxHausdorffDistance MACOSX_BUNDLE milxHausdorffDistance.cpp)
TARGET_LINK_LIBRARIES(milxHausdorffDistance milx-SMILI ${VTK_LIBRARIES} ${ITK_LIBRARIES})

# expand the command line wild-cards in msvc
IF(MSVC)
  SET_TARGET_PROPERTIES(milxHausdorffDistance PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

# Installs the executable
INSTALL(TARGETS milxHausdorffDistance 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime
	)
