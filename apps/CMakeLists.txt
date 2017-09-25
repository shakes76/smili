# CMAKE FILE FOR APPS
cmake_minimum_required(VERSION 3.9)
 
# Recurse into the "apps" subdirectory. This does not actually
# cause another cmake executable to run.
add_subdirectory (milxAssistant)

#################################
# Add necessary sub-directories #
#################################
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
    FIND_PACKAGE(Qt5 COMPONENTS Core Gui Widgets WebEngineWidgets REQUIRED QUIET)
  ELSE(VTK_QT_VERSION VERSION_GREATER "4")
    message("REQUIRES QT5 VERSION 5.9.1 OR GREATER")
  ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

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
