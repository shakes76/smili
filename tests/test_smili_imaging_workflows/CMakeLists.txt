project(test_smili_imaging_workflows)

include_directories(${MILXQT_INCLUDE_DIRS})

add_executable(test_smili_imaging_workflows test_smili_imaging_workflows.cpp)
target_link_libraries(test_smili_imaging_workflows
  milx-Qt
  milx-SMILI
  vtk-ext
  ${QT_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
)

if("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS test_smili_imaging_workflows MODULES ${VTK_LIBRARIES})
endif()

add_test(NAME test_smili_imaging_workflows COMMAND test_smili_imaging_workflows)
