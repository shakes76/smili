PROJECT(prototype_image_planes)

include(find)

INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})

ADD_EXECUTABLE(prototype_image_planes prototype_image_planes.cxx)
TARGET_LINK_LIBRARIES(prototype_image_planes milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})
IF("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS prototype_image_planes MODULES ${VTK_LIBRARIES})
ENDIF()
