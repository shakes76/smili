# CMAKE FILE FOR MILXIMAGEAPP
cmake_minimum_required(VERSION 3.9)

PROJECT(milxImageApp)
cmake_policy(SET CMP0020 NEW)


###############################################
# Find all necessary packages and directories #
###############################################
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})
INCLUDE_DIRECTORIES(${SMILI_INCLUDE_DIRS})
link_directories(${ITK_LIBRARY_DIRS} ${SMILI_BINARY_DIR}/lib)

#~ ADD_DEFINITIONS(-DSMILI_DLL)

#########################################
# Add executable and link the libraries #
#########################################
ADD_EXECUTABLE(milxImageApp MACOSX_BUNDLE milxImageApp.cxx)
TARGET_LINK_LIBRARIES(milxImageApp ${ITK_LIBRARIES} ${ZLIB_LIBRARIES} milx-SMILI)

#######################
# Handle UNIX support #
#######################
IF(UNIX)
  configure_file(milxImageApp.in "milxImageApp" @ONLY)
  
  install(FILES "${PROJECT_BINARY_DIR}/milxImageApp" DESTINATION "../../bin" 
    PERMISSIONS WORLD_READ WORLD_EXECUTE)
ENDIF(UNIX)

##############################################
# Expand the command line wild-cards in MSVC #
##############################################
IF(MSVC)
  SET_TARGET_PROPERTIES(milxImageApp PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

##########################
# Install the executable #
##########################
INSTALL(TARGETS milxImageApp 
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT Runtime
)
