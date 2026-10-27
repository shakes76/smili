PROJECT(imageViewer)

include(find)

INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})
link_directories(${SMILI_BINARY_DIR}/lib)

ADD_EXECUTABLE(imageViewer imageViewer.cxx)
TARGET_LINK_LIBRARIES(imageViewer milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})
IF("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS imageViewer MODULES ${VTK_LIBRARIES})
ENDIF()
