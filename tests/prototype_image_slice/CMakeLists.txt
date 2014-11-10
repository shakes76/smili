cmake_minimum_required(VERSION 2.6)
 
PROJECT(prototype_image_slice)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

SET(SMILI_SOURCE_DIR ../../)

INCLUDE_DIRECTORIES(${SMILI_SOURCE_DIR}/include ${SMILI_SOURCE_DIR}/include/Qt ${SMILI_SOURCE_DIR}/include/vtk-ext ${SMILI_SOURCE_DIR}/include/itkVTKGlue)
link_directories(${SMILI_BINARY_DIR}/lib)
 
ADD_EXECUTABLE(prototype_image_slice prototype_image_slice.cxx)
TARGET_LINK_LIBRARIES(prototype_image_slice milx-SMILI ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${ZLIB_LIBRARIES})
