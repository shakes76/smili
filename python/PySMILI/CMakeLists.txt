#
# This file is part of SMILI.
#
# Add include directories
# Create a library called "milxQt" which includes the source files.
include_directories (${MILXQT_INCLUDE_DIRS} ${ZLIB_INCLUDES} ${SMILI_INCLUDE_DIRS} ${ITK_LIBRARY_DIRS})

# get milx Qt MOCS and Libs for linking later
get_target_property(MILXQT_MOCS milx-Qt MILXQT_MOCS)
get_target_property(MILXQT_INSTALL_BINARY_DIR milx-Qt INTERFACE_LINK_DIRECTORIES)
get_target_property(SMILX_EXE smilx OUTPUT_NAME)

#~ message("MILXQT_MOCS ${MILXQT_MOCS}")

# Need to rebuild the Qt MOCS and UIs etc. for linking
# Add Headers into List
#file(GLOB MILXQT_HEADERS_MOC ${LIBRARY_INCLUDE_PATH}/Qt/*.h)
#file(GLOB MILXQT_SRC_UIS ${SMILI_SOURCE_DIR}/ui/*.ui)
#file(GLOB MILXQT_SRC_UI_HEADERS_MOC ${SMILI_SOURCE_DIR}/ui/*.h)
# build ui_XXX files from the XML-style .ui files
#QT_WRAP_UI(MILXQT_UIS ${MILXQT_SRC_UIS})
# this moc's the above variable and appends to the cxx sources
#QT_WRAP_CPP(MILXQT_MOCS ${MILXQT_HEADERS_MOC})
#QT_WRAP_CPP(MILXQT_UI_MOCS ${MILXQT_SRC_UI_HEADERS_MOC})

if(PYSIDE_MAJOR_VERSION GREATER "2")
    set(PySMILI_pyside6_link_libraries Qt${Qt_VERSION_MAJOR}::OpenGLWidgets)
    set(PySMILI_pyside6_src
        ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkimagedata_wrapper.h
        ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkimagedata_wrapper.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkpoints_wrapper.h
        ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkpoints_wrapper.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkpointset_wrapper.h
        ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkpointset_wrapper.cpp
    )
    set(PySMILI_pyside_xml ${CMAKE_CURRENT_SOURCE_DIR}/typesystem_pyside6_smili.xml)
else()
    set(PySMILI_pyside_xml ${CMAKE_CURRENT_SOURCE_DIR}/typesystem_smili.xml)
endif()

# Auto-Genereate files every class will have his cpp/h files
set(PySMILI_SRC
    # individual classes, must be lower case !!!
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/qvtkopenglnativewidget_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/qvtkopenglnativewidget_wrapper.cpp
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/itk_array_wrapper.h
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/itk_array_wrapper.cpp
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/itk_image_wrapper.h
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/itk_image_wrapper.cpp
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtksmartpointerbase_wrapper.h
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtksmartpointerbase_wrapper.cpp
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtksmartpointer_wrapper.h
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtksmartpointer_wrapper.cpp
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkobjectbase_wrapper.h
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkobjectbase_wrapper.cpp
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkobject_wrapper.h
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkobject_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkdataobject_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkdataobject_wrapper.cpp
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkdataset_wrapper.h
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkdataset_wrapper.cpp
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkdataset_wrapper.h
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkdataset_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkpolydata_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/vtkpolydata_wrapper.cpp
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/qvtkwidget_wrapper.h
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/qvtkwidget_wrapper.cpp
    ${PySMILI_pyside6_src}
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milx_math_wrapper.h
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milx_math_wrapper.cpp
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milx_model_wrapper.h
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milx_model_wrapper.cpp
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milx_file_wrapper.h
    # ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milx_file_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtwindow_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtwindow_wrapper.cpp
#~     ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/labelledaction_wrapper.h
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtunifiedwindow_wrapper.h
    #${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtunifiedwindow_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtrenderwindow_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtrenderwindow_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtimage_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtimage_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtmodel_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtmodel_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtfile_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtfile_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtplot_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtplot_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtmain_wrapper.h
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/milxqtmain_wrapper.cpp
    # global module wrapper
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/pysmili_module_wrapper.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_BINDING_NAMESPACE}/pysmili_python.h

    ${MILXQT_MOCS}
)

# includes necessary to parse and build the classes specified on typesystem
set(PySMILI_include_paths
    #${My_Paths}
    #${MY_CUSTOM_INCLUDE_LIST}
    $<JOIN:$<REMOVE_DUPLICATES:$<TARGET_PROPERTY:milx-Qt,INTERFACE_INCLUDE_DIRECTORIES>>, >
    #$<TARGET_PROPERTY:milx-Qt,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:VTK::CommonCore,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:ITKCommon,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:itkvcl,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:itkvnl_algo,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:itkvnl,INTERFACE_INCLUDE_DIRECTORIES>
    #${ITK_DIR}/../Modules/ThirdParty/VNL/src/vxl/core
    #$<TARGET_PROPERTY:ITKVNLInstantiation,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:VTK::CommonDataModel,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:VTK::RenderingAnnotation,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:VTK::FiltersModeling,INTERFACE_INCLUDE_DIRECTORIES>
    #$<TARGET_PROPERTY:VTK::GUISupportQt,INTERFACE_INCLUDE_DIRECTORIES>
    #${VTK_DIR}/../GUISupport/Qt
    ${MILXQT_INCLUDE_DIRS}
    ${SMILI_INCLUDE_DIRS}
    #${ITK_LIBRARY_DIRS}
)

# A list of paths where shiboken should look for typesystem
set(PySMILI_typesystem_paths
    # PySide path, this variable was exposed by FindPySide2.cmake
    ${PYSIDE_TYPESYSTEMS}
)

# Include flags/path that will be set in 'target_include_directories'
set(PySMILI_target_include_directories ${CMAKE_SOURCE_DIR}/src/Qt)

# Libraries that will be necessary to link the target, this will used in the command 'target_link_libraries'
set(PySMILI_target_link_libraries
    ${VTK_LIBRARIES}
    vtk-ext
    milx-SMILI
    milx-Qt
    Qt${Qt_VERSION_MAJOR}::Core
    Qt${Qt_VERSION_MAJOR}::Gui
    Qt${Qt_VERSION_MAJOR}::Widgets
    ${PySMILI_pyside6_link_libraries}
    # PySide2::pyside2
    ${Python3_LIBRARIES}
)

# changes on these files should trigger a new generation
set(PySMILI_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/smili_global.h
    #${CMAKE_SOURCE_DIR}/include/milxMath.h
    # ${CMAKE_SOURCE_DIR}/include/milxFile.h
    # ${CMAKE_SOURCE_DIR}/include/milxModel.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtWindow.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtRenderWindow.h
    #${CMAKE_SOURCE_DIR}/include/Qt/milxQtUnifiedWindow.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtModel.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtImage.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtFile.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtPlot.h
    ${CMAKE_SOURCE_DIR}/include/Qt/milxQtMain.h
)

create_python_bindings(
  "${PYTHON_BINDING_NAMESPACE}"
  "${PySMILI_typesystem_paths}"
  "${PySMILI_include_paths}"
  "${PySMILI_SRC}"
  "${PySMILI_target_include_directories}"
  "${PySMILI_target_link_libraries}"
  ${CMAKE_CURRENT_SOURCE_DIR}/smili_global.h
  ${PySMILI_pyside_xml}
  "${PySMILI_DEPENDS}"
  ${CMAKE_CURRENT_BINARY_DIR}
)

# Make module import from build dir work
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.cmake ${CMAKE_CURRENT_BINARY_DIR}/__init__.py @ONLY)

#~ get_target_property(QT_CORE_LOCATION Qt6::Core LOCATION)
#~ get_target_property(QT_GUI_LOCATION Qt6::Gui LOCATION)
#~ get_target_property(QT_WIDGETS_LOCATION Qt6::Widgets LOCATION)
#~ get_target_property(QT_WEBENGINEWIDGETS_LOCATION Qt6::WebEngineWidgets LOCATION)
#~ get_target_property(QT_OPENGLWIDGETS_LOCATION Qt6::OpenGLWidgets LOCATION)
#~ get_target_property(QT_OPENGL_LOCATION Qt6::OpenGL LOCATION)
#~ get_target_property(QT_QUICK_LOCATION Qt6::Quick LOCATION)
#~ get_target_property(QT_QUICKWIDGETS_LOCATION Qt6::QuickWidgets LOCATION)
#~ get_target_property(QT_PRINTSUPPORT_LOCATION Qt6::PrintSupport LOCATION)

#~ # Retrieve the INTERFACE_LINK_LIBRARIES property
#~ get_target_property(MILXQT_LINKED_LIBS milx-Qt INTERFACE_LINK_LIBRARIES)
#~ # Print the list of linked libraries
#~ #message(STATUS "milx-Qt Libraries: ${MILXQT_LINKED_LIBS}")

# 1. Initialize an empty list to store the filenames
set(MILXQT_FILENAMES "")
# 2. Iterate through the VTK_LIBRARIES list
foreach(LIB_ITEM ${MILXQT_LINKED_LIBS})
    # Check if the item is a known CMake target (e.g., VTK::CommonCore)
    if(TARGET ${LIB_ITEM})
        # Retrieve the full path to the library file
        get_target_property(LIB_PATH ${LIB_ITEM} IMPORTED_LOCATION_RELEASE)

        if(LIB_PATH)
            # Extract just the filename (e.g., libvtkCommonCore-9.2.so)
            #get_filename_component(LIB_NAME "${LIB_PATH}" NAME)
            list(APPEND MILXQT_FILENAMES "${LIB_PATH}")
        endif()
    else()
        # If it's a raw path instead of a target, extract the name directly
        #get_filename_component(LIB_NAME "${LIB_ITEM}" NAME)
        list(APPEND MILXQT_FILENAMES "${LIB_ITEM}")
    endif()
endforeach()
message("MILXQT FILENAMES: ${MILXQT_FILENAMES}")

# install
message("Install Bindings at ${${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX}")
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
    $<TARGET_FILE:milx-SMILI>
    $<TARGET_FILE:milx-Qt>
#~     $<TARGET_RUNTIME_DLLS:milx-Qt>
    $<TARGET_FILE:vtk-ext>
    ${PYSIDE_LIBRARY}
    ${SHIBOKEN_LIBRARY}
#~     ${QT_CORE_LOCATION}
#~     ${QT_GUI_LOCATION}
#~     ${QT_WIDGETS_LOCATION}
#~     ${QT_WEBENGINEWIDGETS_LOCATION}
#~     ${QT_OPENGLWIDGETS_LOCATION}
#~     ${QT_OPENGL_LOCATION}
#~     ${QT_QUICK_LOCATION}
#~     ${QT_QUICKWIDGETS_LOCATION}
#~     ${QT_PRINTSUPPORT_LOCATION}
#~     ${MILXQT_FILENAMES}
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

##########################
# Install the executable #
##########################
install(TARGETS ${PYTHON_BINDING_NAMESPACE}
        DESTINATION "python"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

INCLUDE(InstallRequiredSystemLibraries)

#set(APPS "${${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX}/${PYTHON_BINDING_NAMESPACE}.pyd")
#~ set(APPS "${CMAKE_CURRENT_BINARY_DIR}/Release/${PYTHON_BINDING_NAMESPACE}.pyd")
#~ set(APPS "${SMILX_EXE}")

# Pull out other dependencies
#SET(LIBRARY_DIRS ${MILXQT_INSTALL_BINARY_DIR};${SMILI_BINARY_DIR}/lib;${VTK_LIBRARY_DIRS};${ITK_LIBRARY_DIRS})

#run qt5 windeploy for Qt platform libs
if(WIN32)
    include(Windeployqt)
    #~ windeployqt(${PYTHON_BINDING_NAMESPACE} ${CMAKE_CURRENT_BINARY_DIR})
    windeployqt(${PYTHON_BINDING_NAMESPACE} ${${PROJECT_NAME}_PYTHON_BINDINGS_INSTALL_PREFIX})
    MESSAGE("Windows Qt${VTK_QT_VERSION} Bindings Deployment")
endif()

#~ install(CODE "
#~ include(BundleUtilities)
#~ fixup_bundle(\"${APPS}\" \"${Plugins}\" \"${MILXQT_INSTALL_BINARY_DIR}\")
#~ " COMPONENT Runtime)
