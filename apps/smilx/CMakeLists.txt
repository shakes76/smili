# CMAKE FILE FOR SMILX
cmake_minimum_required(VERSION 3.9)

#############################
# Handle suport for MAC OSX #
#############################
IF(APPLE) #Assumes MAC OSX
    cmake_minimum_required(VERSION 2.8.8) #for deplotqt module
ELSE(APPLE)
    cmake_minimum_required(VERSION 2.8)
ENDIF(APPLE)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

PROJECT(smilx)
cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0020 NEW)

# Find libraries
include(find)

#######################################################################################
# This sets the windows build that will need the special winmain@16 call. Qt provides #
# this for us in the qtmain.lib file. Using this cmake code will ensure we have it    #
# linked into our build. Not needed on Unix/OS X/Linux which is why we have the       #
# IF(WIN32) conditional. 							      #
#######################################################################################
IF (WIN32)
    SET (QT_USE_QTMAIN TRUE)
ENDIF (WIN32)

###############################
# Find all necessary packages #
###############################
find_package(Qt5 COMPONENTS Widgets REQUIRED QUIET)
set(MILXQT_DIR ${SMILI_SOURCE_DIR})
ADD_DEFINITIONS(-DMILX_DLL) # DLL Export MILX
ADD_DEFINITIONS(-DMILX_PLUGIN_DLL) # DLL Export MILX

#########################
# Add Headers into List #
#########################
file(GLOB MILXQT_HEADERS ${MILXQT_DIR}/include/*.h)
file(GLOB MILXQT_LIBS ${SMILI_BINARY_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB MILXQT_PLUGINS ${SMILI_PLUGINS_PATH}/*${CMAKE_SHARED_LIBRARY_SUFFIX}*)
file(GLOB VTK_QT_LIBS ${VTK_LIBRARY_DIRS}/*Q*${CMAKE_SHARED_LIBRARY_SUFFIX}*)

include_directories (${MILXQT_INCLUDE_DIRS})
link_directories(${MILXQT_INSTALL_BINARY_DIR})

#####################################################################################
# Iconing and APP variable for fixup, needs to be before ADD_EXECUTABLE for Windows #
#####################################################################################
set(APPS "\${CMAKE_INSTALL_PREFIX}/bin/smilx")
IF (WIN32)
    # There is a bug in NSI that does not handle full unix paths properly.
    set(SMILX_ICON_FILE ${PROJECT_SOURCE_DIR}/smilx.rc)
    set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}\\\\smilx.rc")
    set(APPS "\${CMAKE_INSTALL_PREFIX}/bin/smilx.exe")
ELSEIF(APPLE)
    set(APPS "\${CMAKE_INSTALL_PREFIX}/smilx.app")  
    set(SMILX_ICON_FILE ${PROJECT_SOURCE_DIR}/smilx.icns)
    set(MACOSX_BUNDLE_ICON_FILE ${SMILX_ICON_FILE})
    set_source_files_properties(${MACOSX_BUNDLE_ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
    set(SMILX_LD_LIBRARY_PATH "DYLD_LIBRARY_PATH")
ELSE(WIN32)
    set(SMILX_ICON_FILE ${PROJECT_SOURCE_DIR}/smilx.png)
    set(SMILX_LD_LIBRARY_PATH "LD_LIBRARY_PATH")
ENDIF (WIN32)

#########################################
# Add executable and link the libraries #
#########################################
ADD_EXECUTABLE(smilx WIN32 MACOSX_BUNDLE ${SMILX_ICON_FILE} smilx.cpp) #no MACOSX_BUNDLE since custom one created
TARGET_LINK_LIBRARIES(smilx milx-Qt Qt5::Widgets ${VTK_LIBRARIES} ${ITK_LIBRARIES})

##############################################
# Expand the command line wild-cards in MSVC #
##############################################
IF(MSVC)
    SET_TARGET_PROPERTIES(smilx PROPERTIES LINK_FLAGS "setargv.obj")
ENDIF(MSVC)

#############
# Packaging #
#############
SET(CPACK_PACKAGE_NAME "sMILX")
SET(CPACK_PACKAGE_VENDOR "UQ/CSIRO")
SET(CPACK_PACKAGE_VERSION_MAJOR "${SMILX_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${SMILX_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${SMILX_VERSION_PATCH}")
SET(CPACK_MONOLITHIC_INSTALL 1)
SET(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/license.txt)
SET(CPACK_PACKAGE_VERSION ${SMILX_VERSION_STRING})
SET(CPACK_PACKAGE_FILE_NAME "smilx-Setup-v${SMILX_VERSION_STRING}")

##########################
# Install the executable #
##########################
INSTALL(TARGETS smilx 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime
	)

INCLUDE(InstallRequiredSystemLibraries)

IF(BUILD_PLUGINS AND NOT APPLE)
    INSTALL(DIRECTORY ${SMILI_PLUGINS_PATH} DESTINATION . COMPONENT Runtime)
ENDIF(BUILD_PLUGINS AND NOT APPLE)

SET(CPACK_PACKAGE_EXECUTABLES milxOverlay milxOverlay milxImageApp milxImageApp milxModelApp milxModelApp milxDeformableModelApp milxDeformableModelApp milxLabelVisualisation milxLabelVisualisation milxDICOMApp milxDICOMApp smilx sMILX)

IF (WIN32)
    set(_add_executable_param WIN32)
    
    SET(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}\\\\smilx.ico")
    SET(CPACK_NSIS_MUI_UNIICON "${PROJECT_SOURCE_DIR}\\\\smilx.ico")
    # There is a bug in NSI that does not handle full unix paths properly. Make
    # sure there is at least one set of four (4) backlasshes.
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "SMILX")
    SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "smilx ${CPACK_PACKAGE_VERSION}")
    SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\smilx.exe")
    SET(CPACK_NSIS_DISPLAY_NAME "SMILX")
    SET(CPACK_NSIS_HELP_LINK "http://qcplusplus.sourceforge.net/")
    SET(CPACK_NSIS_URL_INFO_ABOUT "http://qcplusplus.sourceforge.net/")
    SET(CPACK_NSIS_CONTACT "shakes76@users.sourceforge.net")
    SET(CPACK_NSIS_MODIFY_PATH ON)
    
    # Pull out other dependencies
    SET(LIBRARY_DIRS ${SMILI_BINARY_DIR}/lib;${MILXQT_INSTALL_BINARY_DIR};${VTK_LIBRARY_DIRS};${ITK_LIBRARY_DIRS};${QT_LIBRARY_DIR};${PYTHONQT_DIR}/lib)    
    
    install(CODE "
    file(GLOB_RECURSE PLUGINS
      \"\${CMAKE_INSTALL_PREFIX}/plugins/*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
    include(BundleUtilities)
    fixup_bundle(\"${APPS}\" \"\${PLUGINS}\" \"${LIBRARY_DIRS}\")
    " COMPONENT Runtime)
    
    # For windows, use nsis
    SET(CPACK_GENERATOR NSIS)
ELSE(WIN32) # Rest of the OSes (Unix and Apple atm)
    SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr/local/sMILX.${CPACK_PACKAGE_VERSION}")
    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "sMILX is a lightweight 3D viewer and graphical frontend for medical image processing.")
    SET(CPACK_PACKAGE_DESCRIPTION "sMILX provides easy ways to do medical image processing and scientific visualisation.")
    SET(CPACK_STRIP_FILES "\${CMAKE_INSTALL_PREFIX}/bin/smilx")
    SET(CPACK_SOURCE_STRIP_FILES "")
  
    configure_file(smilx.sh.in "smilx" @ONLY) #to bin
    configure_file(smilx.sh.in "smilx.sh" @ONLY) #startup
    configure_file(smilx.desktop.in "smilx.desktop" @ONLY)

    IF(APPLE) # Assumes MAC OSX
        SET(CPACK_BUNDLE_NAME "sMILX")
        SET(CPACK_BUNDLE_STARTUP_COMMAND "${PROJECT_BINARY_DIR}/smilx.sh")
        #SET(CPACK_BUNDLE_PLIST "${PROJECT_SOURCE_DIR}/Info-CMake.plist")
        SET(CPACK_BUNDLE_ICON ${SMILX_ICON_FILE})
    
        # Set the OS X Bundle specific CMake variables which will be used to populate the plist for
        # the application bundle
        ConfigureMacOSXBundlePlist(${PROJECT_NAME} "" ${PROJECT_SOURCE_DIR} ${SMILX_VERSION_STRING})

        # Make sure we can find the 'ibtool' program. If we can NOT find it we
        # skip generation of this project
        # find_program(IBTOOL ibtool HINTS "/usr/bin" "${OSX_DEVELOPER_ROOT}/usr/bin")
        # if (${IBTOOL} STREQUAL "IBTOOL-NOTFOUND")
        #   message(SEND_ERROR "ibtool can not be found and is needed to compile the .xib files. It should have been installed with the Apple developer tools. The default system paths were searched in addition to ${OSX_DEVELOPER_ROOT}/usr/bin")
        # endif()

        # Make sure the 'Resources' Directory is correctly created before we build
        # add_custom_command (TARGET smilx PRE_BUILD COMMAND mkdir -p ${EXECUTABLE_OUTPUT_PATH}/\${CONFIGURATION}/smilx.app/Contents/Resources)

        IF(BUILD_PLUGINS)
            message("Plugin Directory is ${SMILI_PLUGINS_PATH}")
            INSTALL(DIRECTORY ${SMILI_PLUGINS_PATH} DESTINATION ${EXECUTABLE_OUTPUT_PATH}/\${CONFIGURATION}/smilx.app/Contents COMPONENT Runtime)
        ENDIF(BUILD_PLUGINS)

        set(CMAKE_INSTALL_PREFIX "")
        SET(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
        set(CPACK_GENERATOR DragNDrop)
        #set(CPACK_GENERATOR DragNDrop; PackageMaker)
        set(CPACK_DMG_FORMAT "UDBZ")
        set(CPACK_DMG_VOLUME_NAME "${CPACK_PACKAGE_NAME}")
        set(CPACK_SYSTEM_NAME "OSX")

        # Pull out other dependencies
        SET(LIBRARY_DIRS ${SMILI_BINARY_DIR}/lib;${MILXQT_INSTALL_BINARY_DIR};${VTK_LIBRARY_DIRS};${ITK_LIBRARY_DIRS};${QT_LIBRARY_DIR};${PYTHONQT_DIR}/lib)    

        IF(CPACK_BINARY_PACKAGEMAKER) #Package Maker 
            set(CMAKE_INSTALL_PREFIX "/Applications") 
            # install(CODE "
            # file(GLOB_RECURSE PLUGINS
            #   \"\${SMILI_PLUGINS_PATH}/*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
            # set(BU_CHMOD_BUNDLE_ITEMS ON)
            # include(BundleUtilities)
            # fixup_bundle(\"${APPS}\" \"\${PLUGINS}\" \"${LIBRARY_DIRS}\")
            # " COMPONENT Runtime)
        # ELSE(CPACK_BINARY_PACKAGEMAKER)
            # Install VTK files 
            # FILE(GLOB vtkFiles "${VTK_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
            # INSTALL(FILES ${vtkFiles} DESTINATION lib COMPONENT Runtime
            # PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

            # Install Qt files
            # FILE(GLOB qtFiles_QtCore "${QT_LIBRARY_DIR}/QtCore.framework/*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
            # INSTALL(FILES ${qtFiles_QtCore} DESTINATION lib COMPONENT Runtime
            # PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)    
            #SET(CPACK_GENERATOR Bundle)
        ENDIF(CPACK_BINARY_PACKAGEMAKER) #Package Maker

        include(DeployQt4)
        set(EXECUTABLE "${PROJECT_NAME}.app")
        FILE(GLOB PLUGINS "${SMILI_PLUGINS_PATH}/*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
        install_qt4_executable("${EXECUTABLE}" "${PLUGINS} qsqlite" "" ${LIBRARY_DIRS})
    ELSE(APPLE) # UNIX etc.
        set(_add_executable_param UNIX)

        # Install all necessary libraries
        install(FILES ${MILXQT_LIBS} DESTINATION lib
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)

        # Install Python QT libraries if built
        IF(BUILD_PYTHON_PLUGIN)
            file(GLOB PYTHONQT_LIBS "${PYTHONQT_DIR}/lib/*PythonQt*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
            install(FILES ${PYTHONQT_LIBS} DESTINATION lib
            PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
        ENDIF(BUILD_PYTHON_PLUGIN)
    
        # Install MSK libraries if built
        IF(BUILD_SSM_PLUGIN)
            FILE(GLOB MSKFiles "${MILXQT_INSTALL_BINARY_DIR}/*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
            install(FILES ${MSKFiles} DESTINATION lib
            PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
        ENDIF(BUILD_SSM_PLUGIN)
    
        # Install MSK libraries if built
        IF(BUILD_FTL_PLUGIN)
            FILE(GLOB FTLFiles1 "${NTTW_DIR}/lib/*frtw*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
            FILE(GLOB FTLFiles2 "${NTTW_DIR}/lib/*nttw*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
            install(FILES ${FTLFiles1} DESTINATION lib
            PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
            install(FILES ${FTLFiles2} DESTINATION lib
            PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
        ENDIF(BUILD_FTL_PLUGIN)
    
        # Install ITK libraries
        # USING PACKAGE MANAGER VERSION - Ubuntu 16.04
#~        IF(ITK_VERSION_MAJOR GREATER 3)
#~            FILE(GLOB itkFiles "${ITK_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
#~        ELSE(ITK_VERSION_MAJOR GREATER 3)
#~            FILE(GLOB itkFiles "${ITK_DIR}/bin/*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
#~        ENDIF(ITK_VERSION_MAJOR GREATER 3)
#~       INSTALL(FILES ${itkFiles} DESTINATION lib COMPONENT Runtime
#~       PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
    
        # Install VTK files 
        IF(VTK_MAJOR_VERSION GREATER 5)
            FILE(GLOB vtkFiles "${VTK_DIR}/lib/*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
        ELSE(VTK_MAJOR_VERSION GREATER 5)
            FILE(GLOB vtkFiles "${VTK_DIR}/bin/*${CMAKE_SHARED_LIBRARY_SUFFIX}*")
        ENDIF(VTK_MAJOR_VERSION GREATER 5)
        INSTALL(FILES ${vtkFiles} DESTINATION lib COMPONENT Runtime
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
    
        # Debian
        SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Shekhar Chandra <shakes76@users.sf.net>")
#~        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libqtcore4, libqtgui4, libqt4-xmlpatterns, libqt4-sql, libqt4-opengl, libqt4-network, libfftw3-3") #10.04
#~        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libqtcore4, libqtgui4, libqt4-xmlpatterns, libqt4-sql, libqt4-opengl, libqt4-network, libvtk5.8, libvtk5.8-qt4, libinsighttoolkit3.20, libpythonqt2.0, libfftw3-3") #12.04
#~        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libqtcore4, libqtgui4, libqt4-xmlpatterns, libqt4-sql, libqt4-opengl, libqt4-network, libvtk5.8, libvtk5.8-qt4, libpythonqt2.0, libfftw3-3") #12.04, no ITK
#~        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libqtcore4, libqtgui4, libqt4-xmlpatterns, libqt4-sql, libqt4-opengl, libqt4-network, libvtk5.8, libvtk5.8-qt4, libpythonqt-dev, libinsighttoolkit4.5, libfftw3-3")
        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libqtcore4, libqtgui4, libqt4-xmlpatterns, libqt4-sql, libqt4-opengl, libqt4-network, libinsighttoolkit4.9, libfftw3-3") #Ubuntu 16.04
#~        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libqtcore4, libqtgui4, libqt4-xmlpatterns, libqt4-sql, libqt4-opengl, libqt4-network, libvtk5.6, libvtk5.6-qt4, libfftw3-3") #11.10

        # RPM
        SET(CPACK_RPM_PACKAGE_SUMMARY ${CPACK_PACKAGE_DESCRIPTION_SUMMARY})
        SET(CPACK_RPM_PACKAGE_NAME ${CPACK_PACKAGE_NAME})
        SET(CPACK_RPM_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION})
        SET(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
        SET(CPACK_RPM_PACKAGE_LICENSE "CSIRO BSD License")
        #SET(CPACK_RPM_PACKAGE_REQUIRES "libqt4, vtk, vtk-qt, libfftw3-3") #opensuse 12.3 and 13.1
        SET(CPACK_RPM_PACKAGE_REQUIRES "qt, vtk, vtk-qt, InsightToolkit, pythonqt, fftw3") #Fedora 21
    
        # Install startup scripts
        install(FILES "${PROJECT_BINARY_DIR}/smilx" DESTINATION "../../bin" 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
        install(FILES "${PROJECT_BINARY_DIR}/smilx.sh" DESTINATION . 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
    
        # Desktop file for menus
        SET(SMILX_DESKTOP_DIR "/usr/share/applications/")
        SET(SMILX_PIXMAPS_DIR "/usr/share/icons/")
#~        FIND_PACKAGE(KDE4)
#~        include(KDE4Internal)
#~        install(FILES ${PROJECT_BINARY_DIR}/smilx.desktop DESTINATION ${XDG_APPS_INSTALL_DIR}) 
        install(FILES ${PROJECT_BINARY_DIR}/smilx.desktop DESTINATION ${SMILX_DESKTOP_DIR}) 
#~        install(FILES ${SMILI_SOURCE_DIR}/resources/smilx_icon.png DESTINATION ${ICON_INSTALL_DIR}) 
        install(FILES ${SMILI_SOURCE_DIR}/resources/smilx_icon.png DESTINATION ${SMILX_PIXMAPS_DIR}) 
#~        kde4_add_app_icon(smilx.cpp "${SMILI_SOURCE_DIR}/resources/smilx_icon.png")
        find_program(XDG-DESKTOP-MENU_EXECUTABLE xdg-desktop-menu)
        INSTALL(CODE "
          execute_process(COMMAND ${XDG-DESKTOP-MENU_EXECUTABLE} install --novendor ${PROJECT_BINARY_DIR}/smilx.desktop)
        ")

        SET(CPACK_GENERATOR TGZ; DEB; RPM)
    ENDIF(APPLE)
  
    set(CPACK_PACKAGE_ICON ${SMILX_ICON_FILE})
ENDIF(WIN32)

INCLUDE(CPack)