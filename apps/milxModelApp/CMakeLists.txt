# CMAKE FILE FOR MILXMODELAPP
cmake_minimum_required(VERSION 3.9)

PROJECT(milxModelApp)
cmake_policy(SET CMP0020 NEW)

###############################################
# Find all necessary packages and directories #
###############################################
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})
link_directories(${VTK_LIBRARY_DIRS} ${SMILI_BINARY_DIR}/lib)

#~ ADD_DEFINITIONS(-DSMILI_DLL)

#########################################
# Add executable and link the libraries #
#########################################
ADD_EXECUTABLE(milxModelApp MACOSX_BUNDLE milxModelApp.cxx)
TARGET_LINK_LIBRARIES(milxModelApp ${VTK_LIBRARIES} milx-SMILI)

#######################
# Handle UNIX support #
#######################
IF(UNIX)
  configure_file(milxModelApp.in "milxModelApp" @ONLY)
  
  install(FILES "${PROJECT_BINARY_DIR}/milxModelApp" DESTINATION "../../bin" 
    PERMISSIONS WORLD_READ WORLD_EXECUTE)
ENDIF(UNIX)

##############################################
# Expand the command line wild-cards in MSVC #
##############################################
IF(MSVC)
  SET_TARGET_PROPERTIES(milxModelApp PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

##########################
# Install the executable #
##########################
INSTALL(TARGETS milxModelApp 
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT Runtime
)
