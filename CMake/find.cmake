#Find the VTK/ITK packages etc.

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

IF(NOT SMILI_FIND_MESSAGE)
    message("Using ITK from ${ITK_DIR}")
ENDIF(NOT SMILI_FIND_MESSAGE)

if(ITK_USE_REVIEW OR "${ITK_VERSION_MAJOR}" GREATER 3)
    add_definitions(-DITK_REVIEW) # used to enable ITK Review features in image class
    IF(NOT SMILI_FIND_MESSAGE)
        message("ITK Review or ITK 4 has been found and will be used.")
    ENDIF(NOT SMILI_FIND_MESSAGE)
    IF( "${ITK_VERSION_MAJOR}" LESS 4 )
        set(ITK_REVIEW_LIBRARIES ITKIOReview)
    ELSE( "${ITK_VERSION_MAJOR}" LESS 4 )
        IF(NOT ITK_USE_REVIEW)
            set(ITK_REVIEW_LIBRARIES )
        ELSE(NOT ITK_USE_REVIEW)
            set(ITK_REVIEW_LIBRARIES ITKReview)
        ENDIF(NOT ITK_USE_REVIEW)
    ENDIF( "${ITK_VERSION_MAJOR}" LESS 4 )
else(ITK_USE_REVIEW OR "${ITK_VERSION_MAJOR}" GREATER 3)
    IF(NOT SMILI_FIND_MESSAGE)
        message("Warning: ITK Review or ITK 4 has not been found. A few image functions may not be available.")
        message("Features missing from imaging: Image Overlay, Overlay Contour")
    ENDIF(NOT SMILI_FIND_MESSAGE)
    set(ITK_REVIEW_LIBRARIES )
endif(ITK_USE_REVIEW OR "${ITK_VERSION_MAJOR}" GREATER 3)

FIND_PACKAGE(VTK)
INCLUDE(${VTK_USE_FILE})

IF(NOT SMILI_FIND_MESSAGE)
    message("Using VTK from ${VTK_DIR}")
ENDIF(NOT SMILI_FIND_MESSAGE)

if(VTK_LIBRARIES)
#~   message("VTK Libraries Defined")
else(VTK_LIBRARIES)
  set(VTK_LIBRARIES vtkHybrid vtkIO)
#~   message("Defined VTK Libraries as ${VTK_LIBRARIES}")
endif(VTK_LIBRARIES)

set(SMILI_FIND_MESSAGE TRUE) #dont keep printing message
