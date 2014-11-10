 
# Recurse into the "apps" subdirectory. This does not actually
# cause another cmake executable to run.
add_subdirectory (milxHausdorffDistanceBetweenSurfaces)
add_subdirectory (milxPlotXY)

if(BUILD_PROTOTYPES)
    add_subdirectory (prototype_plot_xy)
    add_subdirectory (prototype_plot_surface)
    if(VTK_MAJOR_VERSION GREATER 5 OR VTK_MINOR_VERSION GREATER 10)
        add_subdirectory (prototype_image_slice)
    else(VTK_MAJOR_VERSION GREATER 5 OR VTK_MINOR_VERSION GREATER 10)
        message("VTK 5.10.1 or greater not found. Skipping image slice test app.")
    endif(VTK_MAJOR_VERSION GREATER 5 OR VTK_MINOR_VERSION GREATER 10)
    if(ITK_VERSION_MAJOR GREATER 3)
        add_subdirectory (prototype_vectorimage)
    else(ITK_VERSION_MAJOR GREATER 3)
        message("ITK 4 or greater not found. Skipping vector image test app.")
    endif(ITK_VERSION_MAJOR GREATER 3)
endif(BUILD_PROTOTYPES)