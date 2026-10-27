PROJECT(milxHausdorffDistanceBetweenSurfaces)

include(find)

INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})
LINK_DIRECTORIES(${MILXQT_INSTALL_BINARY_DIR})

ADD_EXECUTABLE(milxHausdorffDistanceBetweenSurfaces milxHausdorffDistanceBetweenSurfaces.cxx)
TARGET_LINK_LIBRARIES(milxHausdorffDistanceBetweenSurfaces milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})
IF("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS milxHausdorffDistanceBetweenSurfaces MODULES ${VTK_LIBRARIES})
ENDIF()

# Installs the executable
INSTALL(TARGETS milxHausdorffDistanceBetweenSurfaces 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime
	)
