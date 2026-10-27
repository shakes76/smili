project(test_vtk9_itk5_pipeline)

include_directories(${MILXQT_INCLUDE_DIRS})

add_executable(test_vtk9_itk5_pipeline test_vtk9_itk5_pipeline.cpp)
target_link_libraries(test_vtk9_itk5_pipeline
  milx-Qt
  milx-SMILI
  vtk-ext
  ${QT_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
)

if("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS test_vtk9_itk5_pipeline MODULES ${VTK_LIBRARIES})
endif()

add_test(NAME test_vtk9_itk5_pipeline COMMAND test_vtk9_itk5_pipeline)
