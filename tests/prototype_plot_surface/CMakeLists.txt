PROJECT(prototype_plot_surface)

include(find)

INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})

ADD_EXECUTABLE(prototype_plot_surface prototype_plot_surface.cxx)
TARGET_LINK_LIBRARIES(prototype_plot_surface milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})
IF("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS prototype_plot_surface MODULES ${VTK_LIBRARIES})
ENDIF()
