#
# This file is part of SMILI.
#
# Add include directories
# Create a library called "milxQt" which includes the source files.
include_directories (${MILXQT_INCLUDE_DIRS} ${ZLIB_INCLUDES})

# Auto-Genereate files every class will have his cpp/h files
set(PySMILI_SRC
    # individual classes, must be lower case !!!
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milx_file_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milx_file_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtwindow_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtwindow_wrapper.cpp
#~     ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/labelledaction_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtrenderwindow_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtrenderwindow_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtimage_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtimage_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtmodel_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtmodel_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtfile_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtfile_wrapper.cpp
    # global module wrapper
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/pysmili_module_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/pysmili_python.h
)

# includes necessary to parse and build the classes specified on typesystem
set(PySMILI_include_paths
    $<JOIN:$<TARGET_PROPERTY:milx-Qt,INTERFACE_INCLUDE_DIRECTORIES>,${PATH_SEP}>
)

# A list of paths where shiboken should look for typesystem
set(PySMILI_typesystem_paths
    # PySide path, this variable was exposed by FindPySide2.cmake
    ${PYSIDE_TYPESYSTEMS}
)

# Include flags/path that will be set in 'target_include_directories'
set(PySMILI_target_include_directories
    ${CMAKE_SOURCE_DIR}/src
)

# Libraries that will be necessary to link the target, this will used in the command 'target_link_libraries'
set(PySMILI_target_link_libraries
    vtk-ext
    milx-SMILI
    milx-Qt
    Qt${Qt_VERSION_MAJOR}::Core
    Qt${Qt_VERSION_MAJOR}::Gui
    Qt${Qt_VERSION_MAJOR}::Widgets
    PySide2::pyside2
    ${Python3_LIBRARIES}
)

# changes on these files should trigger a new generation
set(PySMILI_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/smili_global.h
    ${CMAKE_SOURCE_DIR}/include/milxFile.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtWindow.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtRenderWindow.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtModel.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtImage.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtFile.h
)

create_python_bindings(
  "${PYTHON_BINDING_NAMESPACE}"
  "${PySMILI_typesystem_paths}"
  "${PySMILI_include_paths}"
  "${PySMILI_SRC}"
  "${PySMILI_target_include_directories}"
  "${PySMILI_target_link_libraries}"
  ${CMAKE_CURRENT_SOURCE_DIR}/smili_global.h
  ${CMAKE_CURRENT_SOURCE_DIR}/typesystem_smili.xml
  "${PySMILI_DEPENDS}"
  ${CMAKE_CURRENT_BINARY_DIR}
)

# Make module import from build dir work
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.cmake ${CMAKE_CURRENT_BINARY_DIR}/__init__.py @ONLY)

# install
message("Install Bindings at ${${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX}")
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
    $<TARGET_FILE:milx-SMILI>
    $<TARGET_FILE:milx-Qt>
    $<TARGET_FILE:vtk-ext>
    ${PYSIDE_LIBRARY}
    ${SHIBOKEN_LIBRARY}
  DESTINATION
    ${${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX}
)
if(NOT WIN32)
  install(
    FILES
        $<TARGET_LINKER_FILE:milx-SMILI>
        $<TARGET_SONAME_FILE:milx-SMILI>
        $<TARGET_LINKER_FILE:milx-Qt>
        $<TARGET_SONAME_FILE:milx-Qt>
        $<TARGET_LINKER_FILE:vtk-ext>
        $<TARGET_SONAME_FILE:vtk-ext>
        ${PYSIDE_LIBRARY}
        ${SHIBOKEN_LIBRARY}
    DESTINATION
        ${${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX}
  )
endif()
