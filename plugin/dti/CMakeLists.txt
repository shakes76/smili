# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(DTIPlugin CXX)

include(find)

#Find Boost Math libraries, used for sphereical harmonics
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF) 
set(Boost_USE_STATIC_RUNTIME ON)
#~ find_package(Boost REQUIRED COMPONENTS math)
find_package(Boost REQUIRED)

set(SMILI_DIR ${SMILI_SOURCE_DIR})

# Add Headers into List
file(GLOB SMILI_HEADERS ${SMILI_DIR}/include/*.h)

include_directories (${MILXQT_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
link_directories(${SMILI_BINARY_DIR}/lib)

IF(Boost_FOUND)
  ADD_DEFINITIONS("-DHAS_BOOST")
ENDIF(Boost_FOUND)
ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
#ADD_DEFINITIONS(-DQT_NO_DEBUG)
ADD_DEFINITIONS(-DQT_SHARED)
ADD_DEFINITIONS(-DSMILI_DLL)
ADD_DEFINITIONS(-DMILXQT_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_MAKEDLL) #DLL Plugin Export SMILI

# this moc's the above variable and appends to the cxx sources
IF(VTK_QT_VERSION VERSION_GREATER "4")
  QT_WRAP_CPP(SSM_MOCS milxQtDiffusionTensorPlugin.h milxQtDiffusionTensorModel.h milxQtDiffusionTensorImage.h)
ELSE(VTK_QT_VERSION VERSION_GREATER "4")
  QT4_WRAP_CPP(SSM_MOCS milxQtDiffusionTensorPlugin.h milxQtDiffusionTensorModel.h milxQtDiffusionTensorImage.h)
ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

ADD_LIBRARY(DTIPlugin SHARED vtkDiffusionTensorGlyphFilter.cxx milxQtDiffusionTensorPlugin.cpp milxQtDiffusionTensorModel.cpp milxQtDiffusionTensorImage.cpp ${SSM_MOCS})
TARGET_LINK_LIBRARIES(DTIPlugin milx-SMILI milx-Qt ${QT_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${Boost_LIBRARIES})
IF("${VTK_MAJOR_VERSION}" VERSION_GREATER_EQUAL "9")
  vtk_module_autoinit(TARGETS DTIPlugin MODULES ${VTK_LIBRARIES})
ENDIF()

install(TARGETS DTIPlugin
        DESTINATION "plugins"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
