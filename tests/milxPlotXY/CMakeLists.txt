cmake_minimum_required(VERSION 2.6)
 
PROJECT(milxPlotXY)
 
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

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
FIND_PACKAGE( Qt4 REQUIRED )

#--------------------------------------------------------------------------------
# QT_USE_FILE is a variable defined by FIND_PACKAGE( Qt4 ) that contains a path to
# a CMake script that sets up this directory for building with Qt4
INCLUDE( ${QT_USE_FILE} )

# Add Headers into List
file(GLOB MILXQT_HEADERS ${MILXQT_DIR}/include/*.h)
file(GLOB MILXQT_LIBS ${SMILI_BINARY_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)

INCLUDE_DIRECTORIES(${SMILI_SOURCE_DIR}/include ${SMILI_SOURCE_DIR}/include/Qt ${SMILI_SOURCE_DIR}/include/vtk-ext ${SMILI_SOURCE_DIR}/include/itkVTKGlue)
link_directories(${MILXQT_BINARY_DIR})
 
ADD_EXECUTABLE(milxPlotXY milxPlotXY.cxx)
TARGET_LINK_LIBRARIES(milxPlotXY milx-Qt milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})
