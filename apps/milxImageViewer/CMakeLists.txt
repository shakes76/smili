cmake_minimum_required(VERSION 2.6)
 
PROJECT(milxImageViewer)
 
#find libraries
include(find)

#--------------------------------------------------------------------------------
# This sets the windows build that will need the special winmain@16 call. Qt provides
#  this for us in the qtmain.lib file. Using this cmake code will ensure we have it
#  linked into our build. Not needed on Unix/OS X/Linux which is why we have the 
#  IF(WIN32) conditional.
IF (WIN32)
    SET (QT_USE_QTMAIN TRUE)
ENDIF (WIN32)

# Add Headers into List
INCLUDE_DIRECTORIES(${MILXQT_INCLUDE_DIRS})
LINK_DIRECTORIES(${MILXQT_INSTALL_BINARY_DIR})
 
ADD_EXECUTABLE(milxImageViewer milxImageViewer.cpp)
TARGET_LINK_LIBRARIES(milxImageViewer milx-Qt ${QT_LIBRARIES} milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})
