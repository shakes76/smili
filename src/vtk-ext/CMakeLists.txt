# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

PROJECT(vtk-ext)

# Add Sources
#set(VTK_EXT_SOURCES vtkInteractorStyleImage2 vtkImageViewer3 vtkAreaSimplificationMetric)
IF(VTK_MAJOR_VERSION VERSION_GREATER "6") #VTK 7+ deprecated vtkTensor
  set(VTK_EXT_SOURCES vtkInteractorStyleImage2 vtkImageViewer3 vtkTensor)
ELSE(VTK_MAJOR_VERSION VERSION_GREATER "6")
  set(VTK_EXT_SOURCES vtkInteractorStyleImage2 vtkImageViewer3)
ENDIF(VTK_MAJOR_VERSION VERSION_GREATER "6")
set(VTK_EXT_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(VTK_EXT_BINARY_DIR ${PROJECT_BINARY_DIR})
# Add Headers into List
file(GLOB VTK_EXT_HEADERS ${LIBRARY_INCLUDE_PATH}/vtk-ext/*.h)
file(GLOB VTK_EXT_HEADERS_TEMPLATES ${LIBRARY_INCLUDE_PATH}/vtk-ext/*.txx)
file(GLOB VTK_EXT_HEADERS_MOC ${LIBRARY_INCLUDE_PATH}/vtk-ext/*.h)

#Win32 requirement for VTK classes in DLLs
IF(VTK_MAJOR_VERSION VERSION_GREATER "6")
  add_definitions (-DvtkCommonDataModel_EXPORTS)
ENDIF()
IF(VTK_MAJOR_VERSION VERSION_GREATER "5")
  add_definitions (-DvtkInteractionStyle_EXPORTS)
  add_definitions (-DvtkRenderingImage_EXPORTS)
ELSE()
  add_definitions (-DvtkRendering_EXPORTS)
  add_definitions (-DvtkHybrid_EXPORTS)
ENDIF()

include_directories (${LIBRARY_INCLUDE_PATH}/vtk-ext)
link_directories(${VTK_LIBRARY_DIRS} ${LIBRARY_OUTPUT_PATH})

add_library (vtk-ext SHARED ${VTK_EXT_HEADERS} ${VTK_EXT_SOURCES})

# Set the build version (VERSION) and the API version (SOVERSION)
SET_TARGET_PROPERTIES(vtk-ext
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 1)

TARGET_LINK_LIBRARIES(vtk-ext ${VTK_LIBRARIES})

# Installation of the Library
# Installs the header files
install(FILES ${VTK_EXT_HEADERS} ${VTK_EXT_HEADERS_TEMPLATES} DESTINATION "include")
# Installs the library
install(TARGETS vtk-ext
        DESTINATION "lib"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
