# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(animatePlugin CXX)

set(SMILI_DIR ${SMILI_SOURCE_DIR})

# Add Headers into List
file(GLOB SMILI_HEADERS ${SMILI_DIR}/include/*.h)

include_directories (${MILXQT_INCLUDE_DIRS})
link_directories(${SMILI_BINARY_DIR}/lib)

ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
#ADD_DEFINITIONS(-DQT_NO_DEBUG)
ADD_DEFINITIONS(-DQT_SHARED)
ADD_DEFINITIONS(-DMILXQT_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_MAKEDLL) #DLL Plugin Export SMILI

# this moc's the above variable and appends to the cxx sources
IF(VTK_QT_VERSION VERSION_GREATER "4")
    QT5_WRAP_CPP(ANIMATE_MOCS milxQtAnimatePlugin.h milxQtAnimateModel.h)
ELSEIF(VTK_QT_VERSION VERSION_GREATER "4")
    QT4_WRAP_CPP(ANIMATE_MOCS milxQtAnimatePlugin.h milxQtAnimateModel.h)
ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

ADD_LIBRARY(animatePlugin SHARED milxQtAnimateModel.cpp milxQtAnimatePlugin.cpp ${ANIMATE_MOCS})
TARGET_LINK_LIBRARIES(animatePlugin milx-Qt ${QT_LIBRARIES} ${VTK_LIBRARIES})

IF(BUILD_ANIMATE_PLUGIN)
  install(TARGETS animatePlugin
          DESTINATION "plugins"
          PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif(BUILD_ANIMATE_PLUGIN)
