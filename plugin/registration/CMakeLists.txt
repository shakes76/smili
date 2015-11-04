# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
cmake_minimum_required(VERSION 2.6)
PROJECT(registrationPlugin CXX)
set(SMILI_DIR ${SMILI_SOURCE_DIR})

# Add Headers into List
file(GLOB SMILI_HEADERS ${SMILI_DIR}/include/*.h)

# Find the QT
FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtWebkit)
INCLUDE(${QT_USE_FILE})

# Options: library to use
OPTION(USE_NIFTI "Use Nifti registration library" OFF)
OPTION(USE_ELASTIX "Use Elastix registration library" OFF)

# NiftiReg Options
# Nifti reg
IF(USE_NIFTI)
	IF(WIN32)
		IF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) #Win64 environment
			set(NIFTI_REG_DIR $ENV{ProgramW6432}/NiftyReg CACHE STRING "NIFTI REG - Base Directory")
		ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432}) # else use Win32 environment
			set(NIFTI_REG_DIR $ENV{PROGRAMFILES}/NiftyReg CACHE STRING "NIFTI REG - Base Directory")
		ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8 AND EXISTS $ENV{ProgramW6432})
	ENDIF(WIN32)
	IF(UNIX)
		set(NIFTI_REG_DIR /usr/local CACHE PATH "NIFTI REG - Base Directory")
	ENDIF(UNIX)
	MARK_AS_ADVANCED(NIFTI_REG_DIR)
	message("Using NIFTI REG from ${NIFTI_REG_DIR}")
	set(NIFTI_REG_LIBS _reg_aladin _reg_blockMatching _reg_f3d _reg_femTransformation _reg_globalTransformation _reg_KLdivergence _reg_localTransformation _reg_maths _reg_mutualinformation _reg_ReadWriteImage _reg_resampling _reg_ssd _reg_thinPlateSpline _reg_tools reg_nifti reg_nrrd reg_png z png)

	include_directories (${NIFTI_REG_DIR}/include)
	link_directories(${SMILI_BINARY_DIR}/lib ${NIFTI_REG_DIR}/lib)
	
	set(LIBS ${LIBS} ${NIFTI_REG_LIBS})
	add_definitions(-DUSE_NIFTI)
ENDIF(USE_NIFTI)

# Elastix
IF(USE_ELASTIX)
	SET( ELASTIX_BINARY_DIR "" CACHE PATH "Path to elastix binary folder" )
	SET( ELASTIX_USE_FILE ${ELASTIX_BINARY_DIR}/UseElastix.cmake )
	IF( EXISTS ${ELASTIX_USE_FILE} )
	 MESSAGE( STATUS "Including Elastix settings.")
	 INCLUDE( ${ELASTIX_USE_FILE} )
	ENDIF()
	
	file(GLOB ELASTIX_SRC ${SMILI_SOURCE_DIR}/plugin/registration/elastix/*.cxx)
	include_directories(${SMILI_SOURCE_DIR}/plugin/registration/elastix/)
	SET(ELASTIX_LIBS elxCommon elxCore mevisdcmtiff FixedGenericPyramid FixedRecursivePyramid FixedSmoothingPyramid FullSampler GridSampler MultiInputRandomCoordinateSampler RandomSampler RandomCoordinateSampler RandomSamplerSparseMask BSplineInterpolator ReducedDimensionBSplineInterpolator AdvancedMattesMutualInformationMetric AdvancedMeanSquaresMetric AdvancedNormalizedCorrelationMetric TransformBendingEnergyPenalty CorrespondingPointsEuclideanDistanceMetric KNNGraphAlphaMutualInformationMetric NormalizedMutualInformationMetric TransformRigidityPenalty VarianceOverLastDimensionMetric MovingRecursivePyramid MovingSmoothingPyramid AdaptiveStochasticGradientDescent ConjugateGradient FiniteDifferenceGradientDescent FullSearch Powell QuasiNewtonLBFGS RegularStepGradientDescent StandardGradientDescent MultiMetricMultiResolutionRegistration MultiResolutionRegistration MultiResolutionRegistrationWithFeatures BSplineResampleInterpolator ReducedDimensionBSplineResampleInterpolator MyStandardResampler AdvancedAffineTransformElastix AdvancedBSplineTransform DeformationFieldTransform EulerTransformElastix SplineKernelTransform BSplineStackTransform TranslationTransformElastix KNNlib ANNlib)

	set(LIBS ${LIBS} ${ELASTIX_LIBS})
	set(SRC ${SRC} ${ELASTIX_SRC})
	add_definitions(-DUSE_ELASTIX)
ENDIF(USE_ELASTIX)

# ITK
file(GLOB ITK_REG_SRC ${SMILI_SOURCE_DIR}/plugin/registration/itk/*.cxx)
include_directories(${SMILI_SOURCE_DIR}/plugin/registration/itk/)
set(SRC ${SRC} ${ITK_REG_SRC})

# Include MilxQt
include_directories(${MILXQT_INCLUDE_DIRS})

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


set(LIBS milx-SMILI milx-Qt ${QT_LIBRARIES} ${LIBS})
set(SRC ${SSM_MOCS} ${REGISTRATION_UIS} ${REGISTRATION_SRC} ${SRC})

ADD_LIBRARY(registrationPlugin SHARED ${SRC})
include_directories(${ZLIB_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(registrationPlugin ${LIBS})


install(TARGETS registrationPlugin
        DESTINATION "plugins"
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
