 
# Recurse into the "apps" subdirectory. This does not actually
# cause another cmake executable to run.
add_subdirectory (milxAssistant)

IF(USE_VTK)
  add_subdirectory (milxModelApp)
ENDIF(USE_VTK)

IF(USE_ITK)
  add_subdirectory (milxImageApp)
  add_subdirectory (milxDICOMApp)
ENDIF(USE_ITK)

IF(BUILD_MILXQT)
  IF(${ITK_VERSION_MAJOR} VERSION_GREATER 3 OR ITK_USE_REVIEW)
    add_subdirectory (milxLabelVisualisation)
  ELSE(${ITK_VERSION_MAJOR} VERSION_GREATER 3 OR ITK_USE_REVIEW)
    message("Label Visualisation: ITK Version not sufficient to build app. Use ITK 4 or ITK 3 Review")
  ENDIF(${ITK_VERSION_MAJOR} VERSION_GREATER 3 OR ITK_USE_REVIEW)
  add_subdirectory (milxImageViewer)
  add_subdirectory (milxModelViewer)
  add_subdirectory (milxOverlay)
  add_subdirectory (milxBlend)
  if(BUILD_ANIMATE_PLUGIN)
      add_subdirectory(milxAnimate)
  endif(BUILD_ANIMATE_PLUGIN)
ENDIF(BUILD_MILXQT)

IF(USE_ITK AND USE_VTK)
  add_subdirectory (milxHausdorffDistance)
  add_subdirectory (milxDeformableModelApp)
ENDIF(USE_ITK AND USE_VTK)
