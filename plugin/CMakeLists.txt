#Output path
set(LIBRARY_OUTPUT_PATH ${SMILI_PLUGINS_PATH})
#Commodate Apple .app structure
#if(APPLE)
#    set(LIBRARY_OUTPUT_PATH ${SMILI_PLUGINS_PATH})
#endif(APPLE)

if(BUILD_DICOM_PLUGIN)
    add_subdirectory(dicom)
endif(BUILD_DICOM_PLUGIN)

if(BUILD_SSM_PLUGIN)
    add_subdirectory(ssm)
endif(BUILD_SSM_PLUGIN)

if(BUILD_PYTHON_PLUGIN)
    add_subdirectory(python)    
endif(BUILD_PYTHON_PLUGIN)

#~ if(BUILD_DENOISE_PLUGIN)
#~     add_subdirectory(deNoise)
#~ endif(BUILD_DENOISE_PLUGIN)

if(BUILD_ANIMATE_PLUGIN)
    add_subdirectory(animate)
endif(BUILD_ANIMATE_PLUGIN)

if(BUILD_DTI_PLUGIN)
    add_subdirectory(dti)
endif(BUILD_DTI_PLUGIN)

if(BUILD_FTL_PLUGIN)
    add_subdirectory(finitetransform)
endif(BUILD_FTL_PLUGIN)
