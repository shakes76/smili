# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#detect Qt version (depends on VTK-Qt)
#VTK_QT_VERSION is declared in VTKConfig.cmake in VTK 6.1, vtkQt.cmake only available in out-of-source builds and not installed
MESSAGE("VTK-Qt Version: ${VTK_QT_VERSION}")
IF(VTK_QT_VERSION VERSION_GREATER "4")
  MESSAGE("Using Qt 5 for build")
  FIND_PACKAGE(Qt5Widgets REQUIRED QUIET)
  include_directories(${Qt5Widgets_INCLUDE_DIRS})
  add_definitions(${Qt5Widgets_DEFINITIONS})
  qt5_wrap_cpp(QVTKLibMocSrcs ${QVTKMocHeaders})
  set(QT_LIBRARIES ${Qt5Widgets_LIBRARIES})
ELSE(VTK_QT_VERSION VERSION_GREATER "4")
  FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtWebkit)
  INCLUDE(${QT_USE_FILE})
ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

PROJECT(milx-Qt)

# Add Sources
set(MILXQT_SOURCES milxQtConsole milxQtManager milxQtWindow milxQtRenderWindow milxQtImage milxQtModel milxQtPlot milxQtFile milxQtUnifiedWindow milxQtMain)
set(MILXQT_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(MILXQT_BINARY_DIR ${PROJECT_BINARY_DIR})
# Add Headers into List
file(GLOB MILXQT_HEADERS ${LIBRARY_INCLUDE_PATH}/Qt/*.h)
file(GLOB MILXQT_HEADERS_TEMPLATES ${LIBRARY_INCLUDE_PATH}/Qt/*.txx)
file(GLOB MILXQT_HEADERS_MOC ${LIBRARY_INCLUDE_PATH}/Qt/*.h)
file(GLOB MILXQT_SRC_UIS ${SMILI_SOURCE_DIR}/ui/*.ui)
file(GLOB MILXQT_SRC_UI_HEADERS_MOC ${SMILI_SOURCE_DIR}/ui/*.h)
file(GLOB MILXQT_SRC_UI_CPP ${SMILI_SOURCE_DIR}/ui/*.cpp)
file(GLOB MILXQT_SRC_RCCS ${SMILI_SOURCE_DIR}/*.qrc)
#~ message("QRC: ${MILXQT_SRC_RCCS}")

# 64-bit Build?
#set (MILXQT_64 ON) #True
set (MILXQT_64 OFF) #False

#Win32 requirement for VTK classes in DLLs
add_definitions (-DvtkRendering_EXPORTS)

#Add compiler flags
add_definitions (-DMILXQT_DLL -DMILXQT_MAKEDLL)
if (MILXQT_64)
	add_definitions (-DMILXQT_64)
endif (MILXQT_64)

IF(VTK_QT_VERSION VERSION_GREATER "4")
  # generate rules for building source files from the resources
  QT5_ADD_RESOURCES(MILXQT_RCCS ${MILXQT_SRC_RCCS})
  # build ui_XXX files from the XML-style .ui files
  QT5_WRAP_UI(MILXQT_UIS ${MILXQT_SRC_UIS})
  # this moc's the above variable and appends to the cxx sources
  QT5_WRAP_CPP(MILXQT_MOCS ${MILXQT_HEADERS_MOC})
  QT5_WRAP_CPP(MILXQT_UI_MOCS ${MILXQT_SRC_UI_HEADERS_MOC})
ELSE(VTK_QT_VERSION VERSION_GREATER "4")
  # generate rules for building source files from the resources
  QT4_ADD_RESOURCES(MILXQT_RCCS ${MILXQT_SRC_RCCS})
  # build ui_XXX files from the XML-style .ui files
  QT4_WRAP_UI(MILXQT_UIS ${MILXQT_SRC_UIS})
  # this moc's the above variable and appends to the cxx sources
  QT4_WRAP_CPP(MILXQT_MOCS ${MILXQT_HEADERS_MOC})
  QT4_WRAP_CPP(MILXQT_UI_MOCS ${MILXQT_SRC_UI_HEADERS_MOC})
ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

# Create a library called "milxQt" which includes the source files.
# The extension is already found. Any number of sources could be listed here.
include_directories (${MILXQT_INCLUDE_DIRS} ${ZLIB_INCLUDES})
# Add link directories needed to use ITK.
link_directories(${ITK_LIBRARY_DIRS} ${VTK_LIBRARY_DIRS} ${LIBRARY_OUTPUT_PATH})

add_library (milx-Qt SHARED ${MILXQT_HEADERS} ${MILXQT_SOURCES} ${MILXQT_MOCS} ${MILXQT_UI_MOCS} ${MILXQT_UIS} ${MILXQT_SRC_UI_CPP} ${MILXQT_RCCS})

IF(VTK_QT_VERSION VERSION_GREATER "4")
  message("Qt Version: ${Qt5Widgets_VERSION}")
  qt5_use_modules(milx-Qt Core Gui Widgets WebEngineWidgets)
ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

# Set the build version (VERSION) and the API version (SOVERSION)
SET_TARGET_PROPERTIES(milx-Qt
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 1)

#Needed to prevent imporoper linking when VTK 5 is installed side-by-side with VTK 6
IF("${VTK_MAJOR_VERSION}" LESS 6) #If not VTK 6 or above
  set(VTK_QT_LIBRARIES ${VTK_LIBRARIES} QVTK)
ELSE("${VTK_MAJOR_VERSION}" LESS 6) #If VTK 6 or above
  set(VTK_QT_LIBRARIES ${VTK_LIBRARIES})
ENDIF("${VTK_MAJOR_VERSION}" LESS 6) 

TARGET_LINK_LIBRARIES(milx-Qt milx-SMILI vtk-ext ${VTK_LIBRARIES} ${VTK_QT_LIBRARIES} ${QT_LIBRARIES} ${ITK_LIBRARIES} ${ITK_REVIEW_LIBRARIES} ${ZLIB_LIBRARIES})

# Installation of the Library
# Installs the header files
install(FILES ${MILXQT_HEADERS} ${MILXQT_HEADERS_TEMPLATES} DESTINATION "include")
# Installs the library
install(TARGETS milx-Qt
        DESTINATION "lib"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
