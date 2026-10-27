project(test_qt6_smilx_gui)

include_directories(${MILXQT_INCLUDE_DIRS})

add_executable(test_qt6_smilx_gui test_qt6_smilx_gui.cpp)
target_link_libraries(test_qt6_smilx_gui
  milx-Qt
  milx-SMILI
  vtk-ext
  ${QT_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
)

if("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS test_qt6_smilx_gui MODULES ${VTK_LIBRARIES})
endif()

add_test(NAME test_qt6_smilx_gui COMMAND test_qt6_smilx_gui)
