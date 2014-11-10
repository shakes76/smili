#macros to use

# Macro to provide an option only if a set of other variables are ON.
# Example invocation:
#
#  DEPENDENT_OPTION(USE_FOO "Use Foo" ON "USE_BAR;USE_ZOT" OFF)
#
# If both USE_BAR and USE_ZOT are true, this provides an option called
# USE_FOO that defaults to ON.  Otherwise, it sets USE_FOO to OFF.  If
# the status of USE_BAR or USE_ZOT ever changes, any value for the
# USE_FOO option is saved so that when the option is re-enabled it
# retains its old value.
#

MACRO(DEPENDENT_OPTION option doc default depends force)
  IF(${option}_ISSET MATCHES "^${option}_ISSET$")
    SET(${option}_AVAILABLE 1)
    FOREACH(d ${depends})
      IF(NOT ${d})
        SET(${option}_AVAILABLE 0)
      ENDIF(NOT ${d})
    ENDFOREACH(d)

    IF(${option}_AVAILABLE)
      OPTION(${option} "${doc}" "${default}")
      SET(${option} "${${option}}" CACHE BOOL "${doc}" FORCE)
    ELSE(${option}_AVAILABLE)
      IF(NOT ${option} MATCHES "^${option}$")
        SET(${option} "${${option}}" CACHE INTERNAL "${doc}")
      ENDIF(NOT ${option} MATCHES "^${option}$")
      SET(${option} ${force})
    ENDIF(${option}_AVAILABLE)

  ELSE(${option}_ISSET MATCHES "^${option}_ISSET$")
    SET(${option} "${${option}_ISSET}")
  ENDIF(${option}_ISSET MATCHES "^${option}_ISSET$")
ENDMACRO(DEPENDENT_OPTION)

#-------------------------------------------------------------------------------
# This macro will set all the variables necessary to have a "good" OS X Application
# bundle. The variables are as follows:
#  PROJECT_NAME - which can be taken from the ${PROJECT_NAME} variable is needed
#  DEBUG_EXTENSION - The extension used to denote a debug built Application. Typically
#   this is '_debug'
#  ICON_FILE_PATH - The complete path to the bundle icon file
#  VERSION_STRING - The version string that you wish to use for the bundle. For OS X
#   this string is usually XXXX.YY.ZZ in type. Look at the Apple docs for more info
#-------------------------------------------------------------------------------
macro(ConfigureMacOSXBundlePlist PROJECT_NAME DEBUG_EXTENSION ICON_FILE_PATH VERSION_STRING)
    # message(STATUS "ConfigureMacOSXBundlePlist for ${PROJECT_NAME} ")
    IF(CMAKE_BUILD_TYPE MATCHES "Release")
    SET(DBG_EXTENSION "")
    else()
    set(DBG_EXTENSION ${DEBUG_EXTENSION})
    endif()
    get_filename_component(ICON_FILE_NAME "${ICON_FILE_PATH}" NAME)

    #CFBundleGetInfoString
    SET(MACOSX_BUNDLE_INFO_STRING "${PROJECT_NAME}${DBG_EXTENSION} Version ${VERSION_STRING}, Copyright 2014, CSIRO Australia.")
    SET(MACOSX_BUNDLE_ICON_FILE ${ICON_FILE_NAME})
    SET(MACOSX_BUNDLE_GUI_IDENTIFIER "${PROJECT_NAME}${DBG_EXTENSION}")
    #CFBundleLongVersionString
    SET(MACOSX_BUNDLE_LONG_VERSION_STRING "${PROJECT_NAME}${DBG_EXTENSION} Version ${VERSION_STRING}")
    SET(MACOSX_BUNDLE_BUNDLE_NAME ${PROJECT_NAME}${DBG_EXTENSION})
    SET(MACOSX_BUNDLE_SHORT_VERSION_STRING ${VERSION_STRING})
    SET(MACOSX_BUNDLE_BUNDLE_VERSION ${VERSION_STRING})
    SET(MACOSX_BUNDLE_COPYRIGHT "Copyright 2014, CSIRO Australia. All Rights Reserved.")
    # These variables are specific to our plist and are NOT standard CMake variables
    set(MACOSX_BUNDLE_NSMAIN_NIB_FILE "MainMenu")
    set(MACOSX_BUNDLE_NSPRINCIPAL_CLASS "NSApplication")

    SET(${PROJECT_NAME}_PROJECT_SRCS ${${PROJECT_NAME}_PROJECT_SRCS} ${ICON_FILE_PATH})
    SET_SOURCE_FILES_PROPERTIES(${ICON_FILE_PATH} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
endmacro()
