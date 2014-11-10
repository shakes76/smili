cmake_minimum_required(VERSION 2.6)
 
PROJECT(prototype_plot_surface)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

INCLUDE_DIRECTORIES(${SMILI_SOURCE_DIR}/include ${SMILI_SOURCE_DIR}/include/Qt ${SMILI_SOURCE_DIR}/include/vtk-ext ${SMILI_SOURCE_DIR}/include/itkVTKGlue)
 
ADD_EXECUTABLE(prototype_plot_surface prototype_plot_surface.cxx)
TARGET_LINK_LIBRARIES(prototype_plot_surface milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})
