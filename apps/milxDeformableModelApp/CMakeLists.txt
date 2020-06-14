# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(milxDeformableModelApp)

#find libraries
include(find)

INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})

link_directories(${VTK_LIBRARY_DIRS} ${SMILI_BINARY_DIR}/lib)

#~ ADD_DEFINITIONS(-DSMILI_DLL)

ADD_EXECUTABLE(milxDeformableModelApp MACOSX_BUNDLE milxDeformableModelApp.cxx)
TARGET_LINK_LIBRARIES(milxDeformableModelApp ${VTK_LIBRARIES} ${ITK_LIBRARIES} milx-SMILI)

IF(UNIX)
  configure_file(milxDeformableModelApp.in "milxDeformableModelApp" @ONLY)
  
  install(FILES "${PROJECT_BINARY_DIR}/milxDeformableModelApp" DESTINATION "../../bin" 
    PERMISSIONS WORLD_READ WORLD_EXECUTE)
ENDIF(UNIX)

# expand the command line wild-cards in msvc
IF(MSVC)
  SET_TARGET_PROPERTIES(milxDeformableModelApp PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

INSTALL(TARGETS milxDeformableModelApp 
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT Runtime
)
