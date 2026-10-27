PROJECT(prototype_plot_xy)

include(find)

INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})

IF("${VTK_MAJOR_VERSION}" LESS 6)
  set(VTK_CHARTS_LIB vtkCharts)
ELSE()
  set(VTK_CHARTS_LIB "")
ENDIF()

ADD_EXECUTABLE(prototype_plot_xy prototype_plot_xy.cxx)
TARGET_LINK_LIBRARIES(prototype_plot_xy milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES} ${VTK_CHARTS_LIB})
IF("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS prototype_plot_xy MODULES ${VTK_LIBRARIES})
ENDIF()
