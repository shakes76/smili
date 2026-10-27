PROJECT(prototype_image_slice)

include(find)

INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})
link_directories(${SMILI_BINARY_DIR}/lib)

ADD_EXECUTABLE(prototype_image_slice prototype_image_slice.cxx)
TARGET_LINK_LIBRARIES(prototype_image_slice milx-SMILI ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${ZLIB_LIBRARIES})
IF("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS prototype_image_slice MODULES ${VTK_LIBRARIES})
ENDIF()
