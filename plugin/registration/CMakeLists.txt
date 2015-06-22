# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

cmake_minimum_required(VERSION 2.6)

PROJECT(registrationPlugin CXX)

#~ FIND_PACKAGE(VTK REQUIRED)
#~ INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtWebkit)
INCLUDE(${QT_USE_FILE})

set(SMILI_DIR ${SMILI_SOURCE_DIR})

# Add Headers into List
file(GLOB SMILI_HEADERS ${SMILI_DIR}/include/*.h)


#NiftiReg Options		
# Nifti reg
IF(WIN32)
    IF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) #Win64 environment
        set(NIFTI_REG_DIR $ENV{ProgramW6432}/NiftyReg CACHE STRING "NIFTI REG - v1.3.9 - Base Directory")
    ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) # else use Win32 environment
        set(NIFTI_REG_DIR $ENV{PROGRAMFILES}/NiftyReg CACHE STRING "NIFTI REG - v1.3.9 - Base Directory")
    ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432})
ENDIF(WIN32)
IF(UNIX)
    set(NIFTI_REG_DIR /usr/local CACHE PATH "NIFTI REG - v1.3.9 - Base Directory")
ENDIF(UNIX)
MARK_AS_ADVANCED(NIFTI_REG_DIR)
message("Using NIFTI REG 1.3.9 from ${NIFTI_REG_DIR}")
set(NIFTI_REG_LIBS _reg_aladin _reg_blockMatching _reg_f3d _reg_femTransformation _reg_globalTransformation _reg_KLdivergence _reg_localTransformation _reg_maths _reg_mutualinformation _reg_ReadWriteImage _reg_resampling _reg_ssd _reg_thinPlateSpline _reg_tools reg_nifti reg_nrrd reg_png z png)

include_directories (${MILXQT_INCLUDE_DIRS} ${NIFTI_REG_DIR}/include)
link_directories(${SMILI_BINARY_DIR}/lib ${NIFTI_REG_DIR}/lib)


# ADD OPTION FOR EXECUTABLE OR NOT
# ADD_DEFINITIONS(-DQT_SHARED)

# QT DEFINES
ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
ADD_DEFINITIONS(-DQT_SHARED)
ADD_DEFINITIONS(-DMILXQT_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_DLL) #DLL Export SMILI
ADD_DEFINITIONS(-DMILXQT_PLUGIN_MAKEDLL) #DLL Plugin Export SMILI

# UI Building
file(GLOB REGISTRATION_SRC_UIS ${SMILI_SOURCE_DIR}/plugin/registration/*.ui)
file(GLOB REGISTRATION_SRC_HEADERS ${SMILI_SOURCE_DIR}/plugin/registration/*.h)
file(GLOB REGISTRATION_SRC ${SMILI_SOURCE_DIR}/plugin/registration/*.cpp)
QT4_WRAP_UI(REGISTRATION_UIS ${REGISTRATION_SRC_UIS})
QT4_WRAP_CPP(SSM_MOCS ${REGISTRATION_SRC_HEADERS})


 
ADD_LIBRARY(registrationPlugin SHARED ${SSM_MOCS} ${REGISTRATION_UIS} ${REGISTRATION_SRC})
include_directories(${ZLIB_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(registrationPlugin milx-SMILI milx-Qt ${QT_LIBRARIES} ${NIFTI_REG_LIBS})


install(TARGETS registrationPlugin
        DESTINATION "plugins"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
