# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(milx-SMILI)

#load modules
#-----------------------------------------------------------------------------
include(zlib)
include(find)

# Add module library
set(SMILI_INC ${LIBRARY_INCLUDE_PATH})

# Add Sources
IF(USE_ITK)
set(SMILI_SOURCES
    milxFile.cxx
   )
ENDIF(USE_ITK)
IF(USE_VTK)
set(SMILI_SOURCES
    ${SMILI_SOURCES}
    utils/milxColourMap.cxx
    milxModel.cxx
   )
ENDIF(USE_VTK)
IF(USE_VTK AND USE_ITK)
set(SMILI_SOURCES
    ${SMILI_SOURCES}
    milxDeformableModel.cxx
   )
ENDIF(USE_VTK AND USE_ITK)

# Add Headers into List
file(GLOB SMILI_HEADERS ${LIBRARY_INCLUDE_PATH}/*.h)
IF(USE_ITK)
  set(SMILI_HEADERS
    ${SMILI_HEADERS}
    #itk extensions
    ${SMILI_INC}/itk-ext/itkVectorShiftScaleImageFilter.h
    ${SMILI_INC}/itk-ext/itkVectorShiftScaleImageFilter.hxx
    )
ENDIF(USE_ITK)
IF(USE_VTK)
  set(SMILI_HEADERS
    ${SMILI_HEADERS}
    ${SMILI_INC}/utils/milxLUT.h
    ${SMILI_INC}/utils/milxColourMap.h
    #itkVTKGlue
    ${SMILI_INC}/itkVTKGlue/itkImageToVTKImageFilter.h
    ${SMILI_INC}/itkVTKGlue/itkImageToVTKImageFilter.txx
    ${SMILI_INC}/itkVTKGlue/itkMeshToVTKPolyData.h
    ${SMILI_INC}/itkVTKGlue/itkMeshToVTKPolyData.txx
    ${SMILI_INC}/itkVTKGlue/itkVTKImageToImageFilter.h
    ${SMILI_INC}/itkVTKGlue/itkVTKImageToImageFilter.txx
    ${SMILI_INC}/itkVTKGlue/itkVTKPolyDataToMesh.h
    ${SMILI_INC}/itkVTKGlue/itkVTKPolyDataToMesh.txx
    )
ENDIF(USE_VTK)

# 64-bit Build?
#set (SMILI_64 ON) #True
set (SMILI_64 OFF) #False

#Add compiler flags
#add_definitions (-DSMILI_DLL -DSMILI_MAKEDLL)
if (SMILI_64)
	add_definitions (-DSMILI_64)
endif (SMILI_64)

set(BUILD_TYPE STATIC)
if(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DSMILI_DLL)
  ADD_DEFINITIONS(-DSMILI_MAKEDLL)
  #Win32 requirement for VTK classes in DLLs
  add_definitions (-DvtkCommon_EXPORTS)
  set(BUILD_TYPE SHARED)
endif(BUILD_SHARED_LIBS)

# Add include directories
include_directories (${SMILI_INC} ${LIBRARY_INCLUDE_PATH} ${SMILI_INCLUDE_DIRS} ${ZLIB_INCLUDES})
# Add link directories needed to use ITK.
link_directories(${ITK_LIBRARY_DIRS} ${VTK_LIBRARY_DIRS})

add_library(milx-SMILI ${BUILD_TYPE} ${SMILI_HEADERS} ${SMILI_SOURCES})

# Set the build version (VERSION) and the API version (SOVERSION)
SET_TARGET_PROPERTIES(milx-SMILI
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 1)

TARGET_LINK_LIBRARIES(milx-SMILI vtk-ext ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${ZLIB_LIBRARIES})

# Installation of the Library
# Installs the header files
install(FILES ${SMILI_HEADERS} DESTINATION "include")
# Installs the library
install(TARGETS milx-SMILI
        DESTINATION "lib"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
        
#build the MILX Qt library if requested
IF(BUILD_VTK_EXTENSIONS AND USE_VTK)
    ADD_SUBDIRECTORY(vtk-ext)
ENDIF(BUILD_VTK_EXTENSIONS AND USE_VTK)

IF(BUILD_MILXQT)
    SET(USE_ITK ON)
    SET(USE_VTK ON)
    ADD_SUBDIRECTORY(Qt)
ENDIF(BUILD_MILXQT)
