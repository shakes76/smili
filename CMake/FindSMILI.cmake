 
# - Try to find SMILI
# Once done this will define
#  SMILI_FOUND - System has SMILI
#  SMILI_INCLUDE_DIRS - The SMILI include directories
#  SMILI_LIBRARIES - The libraries needed to use SMILI
#  SMILI_DEFINITIONS - Compiler switches required for using SMILI

#~ find_package(PkgConfig)
#~ pkg_check_modules(PC_LIBXML QUIET libmilx-SMILI)
#~ set(SMILI_DEFINITIONS ${PC_LIBXML_CFLAGS_OTHER})

find_path(SMILI_INCLUDE_DIR NAMES milxImage.h milxFile.h
          HINTS "~/Dev/Install/include" "/usr/local/sMILX.1.0.0/include" "/usr/local/sMILX.0.99.9/include" "C:/Program Files (x86)/SMILX/include" "C:/Program Files/SMILX/include"
          PATH_SUFFIXES smili )
          
find_path(MILXQT_INCLUDE_DIRS NAMES milxQtImage.h milxQtMain.h
          HINTS "~/Dev/Install/include" "/usr/local/sMILX.1.0.0/include" "/usr/local/sMILX.0.99.9/include" "C:/Program Files (x86)/SMILX/include" "C:/Program Files/SMILX/include"
          PATH_SUFFIXES smili Qt )

find_library(SMILI_LIBRARY NAMES milx-SMILI
             HINTS "~/Dev/Install/lib" "/usr/local/sMILX.1.0.0/lib" "/usr/local/sMILX.0.99.9/lib" "C:/Program Files (x86)/SMILX/lib" "C:/Program Files/SMILX/lib")
find_library(MILXQT_LIBRARY NAMES milx-Qt
             HINTS "~/Dev/Install/lib" "/usr/local/sMILX.1.0.0/lib" "/usr/local/sMILX.0.99.9/lib" "C:/Program Files (x86)/SMILX/lib" "C:/Program Files/SMILX/lib")
find_library(VTK_EXT_LIBRARY NAMES vtk-ext
             HINTS "~/Dev/Install/lib" "/usr/local/sMILX.1.0.0/lib" "/usr/local/sMILX.0.99.9/lib" "C:/Program Files (x86)/SMILX/lib" "C:/Program Files/SMILX/lib")

set(SMILI_LIBRARIES ${SMILI_LIBRARY} ${MILXQT_LIBRARY} ${VTK_EXT_LIBRARY})
MESSAGE("SMILI LIBRARIES: ${SMILI_LIBRARIES}")
set(SMILI_INCLUDE_DIRS ${SMILI_INCLUDE_DIR} ${MILXQT_INCLUDE_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SMILI_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SMILI  DEFAULT_MSG
                                  SMILI_LIBRARY MILXQT_LIBRARY VTK_EXT_LIBRARY SMILI_INCLUDE_DIR)

mark_as_advanced(SMILI_INCLUDE_DIR SMILI_LIBRARIES )
