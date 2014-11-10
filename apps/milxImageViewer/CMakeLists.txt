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

#--------------------------------------------------------------------------------
# this command finds Qt4 libraries and sets all required variables
# note that it's Qt4, not QT4 or qt4
FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtWebkit)

#--------------------------------------------------------------------------------
# QT_USE_FILE is a variable defined by FIND_PACKAGE( Qt4 ) that contains a path to
# a CMake script that sets up this directory for building with Qt4
INCLUDE( ${QT_USE_FILE} )

# Add Headers into List
INCLUDE_DIRECTORIES(${MILXQT_INCLUDE_DIRS})
LINK_DIRECTORIES(${MILXQT_INSTALL_BINARY_DIR})
 
ADD_EXECUTABLE(milxImageViewer milxImageViewer.cpp)
TARGET_LINK_LIBRARIES(milxImageViewer milx-Qt ${QT_LIBRARIES} milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})
