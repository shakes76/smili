cmake_minimum_required(VERSION 2.6)
 
PROJECT(milxHausdorffDistanceBetweenSurfaces)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})
LINK_DIRECTORIES(${MILXQT_INSTALL_BINARY_DIR})
 
ADD_EXECUTABLE(milxHausdorffDistanceBetweenSurfaces milxHausdorffDistanceBetweenSurfaces.cxx)
TARGET_LINK_LIBRARIES(milxHausdorffDistanceBetweenSurfaces milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})

# Installs the executable
INSTALL(TARGETS milxHausdorffDistanceBetweenSurfaces 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime
	)