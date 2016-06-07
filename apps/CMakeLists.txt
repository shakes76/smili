 
# Recurse into the "apps" subdirectory. This does not actually
# cause another cmake executable to run.

IF(USE_VTK)
  add_subdirectory (milxModelApp)
ENDIF(USE_VTK)

IF(USE_ITK)
  add_subdirectory (milxImageApp)
  add_subdirectory (milxDICOMApp)
ENDIF(USE_ITK)

IF(BUILD_MILXQT)
  set(CMAKE_AUTOMOC ON)
  #--------------------------------------------------------------------------------
  # this command finds Qt4 libraries and sets all required variables
  IF(VTK_QT_VERSION VERSION_GREATER "4")
    FIND_PACKAGE(Qt5 COMPONENTS Core Gui Widgets WebView REQUIRED QUIET)
  ELSE(VTK_QT_VERSION VERSION_GREATER "4")
    FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtWebkit)
    # QT_USE_FILE is a variable defined by FIND_PACKAGE( Qt4 ) that contains a path to
    # a CMake script that sets up this directory for building with Qt4
    INCLUDE( ${QT_USE_FILE} )
  ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

  add_subdirectory (milxAssistant)

  IF(${ITK_VERSION_MAJOR} VERSION_GREATER 3 OR ITK_USE_REVIEW)
    add_subdirectory (milxLabelVisualisation)
  ELSE(${ITK_VERSION_MAJOR} VERSION_GREATER 3 OR ITK_USE_REVIEW)
    message("Label Visualisation: ITK Version not sufficient to build app. Use ITK 4 or ITK 3 Review")
  ENDIF(${ITK_VERSION_MAJOR} VERSION_GREATER 3 OR ITK_USE_REVIEW)
  add_subdirectory (milxImageViewer)
  add_subdirectory (milxModelViewer)
  add_subdirectory (milxOverlay)
  if(BUILD_ANIMATE_PLUGIN)
      add_subdirectory(milxAnimate)
  endif(BUILD_ANIMATE_PLUGIN)
ENDIF(BUILD_MILXQT)

IF(USE_ITK AND USE_VTK)
  add_subdirectory (milxHausdorffDistance)
  add_subdirectory (milxDeformableModelApp)
ENDIF(USE_ITK AND USE_VTK)
