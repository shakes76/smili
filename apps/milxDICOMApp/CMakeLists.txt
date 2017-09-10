# CMAKE FILE FOR MILXDICOMAPP
cmake_minimum_required(VERSION 3.9)

PROJECT(milxDICOMApp)
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
ADD_EXECUTABLE(milxDICOMApp MACOSX_BUNDLE milxDICOMApp.cxx)
TARGET_LINK_LIBRARIES(milxDICOMApp ${ITK_LIBRARIES} ${ZLIB_LIBRARIES} milx-SMILI)

#######################
# Handle UNIX support #
#######################
IF(UNIX)
  configure_file(milxDICOMApp.in "milxDICOMApp" @ONLY)
  
  install(FILES "${PROJECT_BINARY_DIR}/milxDICOMApp" DESTINATION "../../bin" 
    PERMISSIONS WORLD_READ WORLD_EXECUTE)
ENDIF(UNIX)

##############################################
# Expand the command line wild-cards in MSVC #
##############################################
IF(MSVC)
  SET_TARGET_PROPERTIES(milxDICOMApp PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

##########################
# Install the executable #
##########################
INSTALL(TARGETS milxDICOMApp 
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT Runtime
)
