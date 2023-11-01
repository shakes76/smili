cmake_minimum_required(VERSION 3.2)
 
PROJECT(SMILI)

#Set output paths
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(LIBRARY_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/include)
set(SMILI_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(SMILI_BINARY_DIR ${PROJECT_BINARY_DIR})
set(CMAKE_MODULE_PATH "${SMILI_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

#load modules
#-----------------------------------------------------------------------------
include(zlib)
include(find)
include(macros)
include(version)

#Add compiler flags
#GCC
if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fwrapv")
endif(CMAKE_COMPILER_IS_GNUCC)
#MSVC
if(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj")
endif(MSVC)

# Enable/Disable Build the testing
OPTION(BUILD_SHARED_LIBS "Build the Libraries as shared libraries" ON)
OPTION(BUILD_GUI_RELATED "Build the MILX Qt Related Components" ON)
OPTION(BUILD_VIEWER "Build the MILX Qt Viewer Application" ON)
OPTION(BUILD_APPLICATIONS "Build the applications for SMILI. Use BUILD_GUI_RELATED to enable Qt applications." ON)
OPTION(BUILD_TESTS "Build the tests for SMILI" OFF)
OPTION(BUILD_PROTOTYPES "Build the prototype tests for SMILI" OFF)

#-----------------------------------------------------------------------------
# Enable/Disable Freesurfer support
OPTION(BUILD_FREESURFER_SUPPORT "Build the Freesurfer support parts of SMILI." OFF)

#dependent options, triggered if another option is set
DEPENDENT_OPTION(USE_ITK "Build all the ITK components" ON "ITK_LIBRARIES" OFF)
DEPENDENT_OPTION(USE_VTK "Build all the VTK components" ON "VTK_LIBRARIES" OFF)
DEPENDENT_OPTION(BUILD_VTK_EXTENSIONS "Build the VTK Extensions Library" ON "USE_VTK" OFF)
DEPENDENT_OPTION(BUILD_MILXQT "Build the MILX Qt Library" ON "BUILD_GUI_RELATED" OFF)
DEPENDENT_OPTION(BUILD_PLUGINS "Build the plugins for the Viewer" OFF "BUILD_VIEWER" OFF)
OPTION(PYTHON_BINDINGS "Build python bindings" OFF)

IF(BUILD_MILXQT)
    #detect Qt version (depends on VTK-Qt)
    #VTK_QT_VERSION is declared in VTKConfig.cmake in VTK 6.1, vtkQt.cmake only available in out-of-source builds and not installed
    MESSAGE("VTK-Qt Version: ${VTK_QT_VERSION}")
    IF(VTK_QT_VERSION GREATER "4")
        find_package("Qt${VTK_QT_VERSION}" COMPONENTS Core Gui Widgets WebEngineWidgets REQUIRED)
#~         qt_standard_project_setup()
        include_directories(${Qt${QT_VERSION_MAJOR}Core_INCLUDE_DIRS} ${Qt${QT_VERSION_MAJOR}WebEngineWidgets_INCLUDE_DIRS})
        add_definitions(${Qt${QT_VERSION_MAJOR}Core_DEFINITIONS})
        qt_wrap_cpp(QVTKLibMocSrcs ${QVTKMocHeaders})
        set(QT_LIBRARIES Qt${QT_VERSION_MAJOR}::Core Qt${QT_VERSION_MAJOR}::WebEngineWidgets)
        set(${PROJECT_NAME}_QT${VTK_QT_VERSION} TRUE)
    ELSE(VTK_QT_VERSION GREATER "4")
        FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtWebkit)
        INCLUDE(${QT_USE_FILE})
    ENDIF(VTK_QT_VERSION GREATER "4")
ENDIF(BUILD_MILXQT)
 
IF(BUILD_PLUGINS)
    OPTION(BUILD_PYTHON_PLUGIN "Build the Python plugin for MILX Qt Viewer" OFF)
    OPTION(BUILD_DICOM_PLUGIN "Build the DICOM plugin for MILX Qt Viewer" OFF)
    OPTION(BUILD_SSM_PLUGIN "Build the Shape Model plugin for MILX Qt Viewer" OFF)
    OPTION(BUILD_DENOISE_PLUGIN "Build the Denoise plugin for MILX Qt Viewer" OFF)
    OPTION(BUILD_ANIMATE_PLUGIN "Build the Animate plugin for MILX Qt Viewer" OFF)
    OPTION(BUILD_DTI_PLUGIN "Build the Diffusion Tensor plugin for MILX Qt Viewer" OFF)
    OPTION(BUILD_FTL_PLUGIN "Build the Finite Transform Library plugin for MILX Qt Viewer" OFF)
    OPTION(BUILD_REGISTRATION_PLUGIN "Build the Image Registration Plugin for MILX Qt " OFF)
ENDIF(BUILD_PLUGINS)

IF(BUILD_FREESURFER_SUPPORT)
    add_definitions(-DFREESURFER_SUPPORT)  
ENDIF(BUILD_FREESURFER_SUPPORT)

#Check ITK/VTK statuses
IF(USE_ITK AND NOT USE_VTK)
  add_definitions(-DITK_ONLY)
  message("Disabling non-ITK components. Enable USE_VTK to enable VTK components.")
ENDIF(USE_ITK AND NOT USE_VTK)
IF(USE_VTK AND NOT USE_ITK)
  add_definitions(-DVTK_ONLY)
  message("Disabling non-VTK components. Enable USE_ITK to enable ITK components.")
ENDIF(USE_VTK AND NOT USE_ITK)
IF(USE_VTK)
  message("VTK Version: ${VTK_VERSION}")
  add_definitions( -DVTK_MAJOR_VERSION=${VTK_MAJOR_VERSION})
ENDIF(USE_VTK)
IF(USE_ITK)
  message("ITK Version: ${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}")
ENDIF(USE_ITK)

#advanced options
MARK_AS_ADVANCED(BUILD_PROTOTYPES)
MARK_AS_ADVANCED(BUILD_VTK_EXTENSIONS)

#Important variables for configurations
set(SMILI_INCLUDE_DIRS ${LIBRARY_INCLUDE_PATH} ${LIBRARY_INCLUDE_PATH}/vtk-ext ${LIBRARY_INCLUDE_PATH}/itk-ext ${LIBRARY_INCLUDE_PATH}/itkVTKGlue ${LIBRARY_INCLUDE_PATH}/utils)
set(MILXQT_INCLUDE_DIRS ${SMILI_INCLUDE_DIRS} ${LIBRARY_INCLUDE_PATH}/Qt ${SMILI_SOURCE_DIR}/ui)
set(SMILI_PLUGINS_PATH ${SMILI_BINARY_DIR}/plugins)

IF(APPLE)
    #Release RPATH handling for MAC OSX
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_MACOSX_RPATH OFF)
    MESSAGE("MAC OSX RPATH STUFF")
    cmake_policy(SET CMP0042 NEW)
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
    #SET(CMAKE_SKIP_RPATH TRUE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # the RPATH to be used when installing, but only if it's not a system directory
    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
       SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    ENDIF("${isSystemDir}" STREQUAL "-1")
ENDIF(APPLE)

# Recurse into the "apps" subdirectory. This does not actually
# cause another cmake executable to run.
add_subdirectory (src)

IF(PYTHON_BINDINGS)
  IF(CMAKE_BUILD_TYPE MATCHES "^[Dd]eb" OR ${PROJECT_NAME}_STATIC)
    MESSAGE(FATAL_ERROR "** Python Bindings are disabled in debug or static builds.")
  ENDIF()
ENDIF()
IF(PYTHON_BINDINGS)
  add_subdirectory (python)
ENDIF(PYTHON_BINDINGS)

IF(BUILD_VIEWER)
    add_subdirectory (apps/smilx)
ENDIF(BUILD_VIEWER)

IF(BUILD_PLUGINS)
    add_subdirectory (plugin)
ENDIF(BUILD_PLUGINS)

IF(BUILD_APPLICATIONS)
    add_subdirectory (apps)
ENDIF(BUILD_APPLICATIONS)

IF(BUILD_TESTS)
    add_subdirectory (tests)
ENDIF(BUILD_TESTS)
