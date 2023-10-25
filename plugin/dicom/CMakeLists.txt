###########################
# DICOM PLUGIN CMAKE FILE #
###########################
PROJECT(DICOMPlugin CXX)
cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0020 NEW)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

##################
# Find libraries #
##################
include(find)

if("${ITK_VERSION_MAJOR}" LESS 4)
    message("DICOM Plugin requires ITK 4 or above for DICOM RT support. Disabling structure set support.")
endif("${ITK_VERSION_MAJOR}" LESS 4)

set(SMILI_DIR ${SMILI_SOURCE_DIR})

#########################
# Add Headers into List #
#########################
file(GLOB SMILI_HEADERS ${SMILI_DIR}/include/*.h)

################################
# Include and link directories #
################################
include_directories (${MILXQT_INCLUDE_DIRS})
link_directories(${SMILI_BINARY_DIR}/lib)

ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
#ADD_DEFINITIONS(-DQT_NO_DEBUG)
ADD_DEFINITIONS(-DQT_SHARED)
ADD_DEFINITIONS(-DMILXQT_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_MAKEDLL) #DLL Plugin Export SMILI
#Win32 requirement for VTK classes in DLLs
#~ ADD_DEFINITIONS (-DvtkFiltering_EXPORTS)

# this moc's the above variable and appends to the cxx sources
IF(VTK_QT_VERSION VERSION_GREATER "4")
  QT_WRAP_CPP(MILXQT_MOCS ${MILXQT_HEADERS_MOC})
ELSE(VTK_QT_VERSION VERSION_GREATER "4")
  QT4_WRAP_CPP(SSM_MOCS milxQtDICOMPlugin.h)
ENDIF(VTK_QT_VERSION VERSION_GREATER "4")

###################################
# Add and link external libraries #
###################################
ADD_LIBRARY(DICOMPlugin SHARED milxQtDICOMPlugin.cpp ${SSM_MOCS})
TARGET_LINK_LIBRARIES(DICOMPlugin milx-SMILI milx-Qt ${QT_LIBRARIES})

##########################
# Install the executable #
##########################
if(NOT APPLE)
	install(TARGETS DICOMPlugin
	        DESTINATION "plugins"
	        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
endif()