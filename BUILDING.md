![SMILI Logo](browse/resources/smili_logo.png?raw=true)

# Building SMILI
* You will need CMake to build SMILI installed on your OS. Get it from cmake.org.
* You will need the Insight Toolkit (ITK) built and/or installed on your OS. ITK 4.5.x or above for the DICOM plugin and vector image support. ITK 3 is reccomended for the Shape Model (SSM) plugin
* You will need Qt 4 Framework. Qt 5 is not yet supported.
* You will need the Visualisation Toolkit (VTK) built and/or installed on your OS. VTK 6.1 or above is reccomended and VTK 5.8 for the SSM plugin.

## Linux
Just install CMake, ITK 4.x and VTK 6.x from your distro's repository. The remaining aspects are covered in a video on youtube found on the SMILI channel:
https://www.youtube.com/channel/UCD-hU6IF2qGlz7roexAUj1Q

### Building SMILI
* From the SMILI source directory create a build directory, change to it and run CMake Curses GUI

  '''
  mkdir build
  cd build
  ccmake ..
  '''
* Select your build environment/compiler with the -G option when you run 'ccmake ..' if you want to use other IDEs like Codeblocks etc.
* Configure and Generate as no other specific settings are required.
* Open solution etc. and build making sure the architecture and build type (Release etc.) is consistent throughout.

#### Python Plugin
This plugin requires PythonQt 2.1 or above to be installed from the repository.

## Windows

### Get Qt 4 Setup
It is best to grab the Qt 4.6 or above binary from the Qt-x64 project of SourceForge (http://sourceforge.net/projects/qtx64/).

### Building ITK
Consult the ITK documentation: http://www.itk.org/Wiki/ITK/Configuring_and_Building
* Download ITK from itk.org. 
* Open CMake-GUI and select the source firectory and create a directory for the binaries
* Select your build environment/compiler
* Turn BUILD_SHARED_LIBS to ON and testing and examples OFF.
* Configure and Generate as no other specific settings are required.
* Open solution etc. and build making sure the architecture and build type (Release etc.) is consistent throughout.

### Building VTK
Consult the VTK documentation: http://www.vtk.org/Wiki/VTK/Configure_and_Build
* Download VTK from vtk.org
* Open CMake-GUI and select the source firectory and create a directory for the binaries
* Select your build environment/compiler
* Turn BUILD_SHARED_LIBS to ON and testing and examples OFF.
* Turn VTK_Group_Qt to ON. It may not find Qt even if you've got it from the Qt-x64 project, just point the QT_QMAKE_EXECUTABLE to the qmake.exe file in the unzipped binary bin directory.
* If wanting to build the Animate plugin, you will need to turn Module_vtkFFMPEG to ON. You will need the FFMPEG binaries. Get the shared and dev files from http://ffmpeg.zeranoe.com/builds/ and fill in the required fields when VTK cant find them.
   You may get issues with the FFMPEG API changing every release. A nice summary of changes can be found at http://sgros.blogspot.com.au/2013/01/deprecated-functions-in-ffmpeg-library.html
* Configure and Generate as no other specific settings are required.
* Open solution etc. and build making sure the architecture and build type (Release etc.) is consistent throughout.

### Building SMILI
* Open CMake-GUI and select the source firectory and create a directory for the binaries
* Select your build environment/compiler
* Qt 4 library will likely not be found so just point the QT_QMAKE_EXECUTABLE to the qmake.exe file in the unzipped binary bin directory.
   If you're only trying to build the GUI-independent sub-library milxSMILI, then just set BUILD_VIEWER to OFF.
* CMake should automatically find the previous VTK and ITK libraries you built. If not, change the ITK_DIR and VTK_DIR appropriately.
* Set the ZLIB variable manually on Windows, just reuse the ITK ZLIB libraries and headers. Otherwise you'll get the error "zlib.h" not found.
* Configure and Generate as no other specific settings are required.
* Open solution etc. and build making sure the architecture and build type (Release etc.) is consistent throughout.

###Building Plugins
#### DICOM Plugin
This plugin just requires ITK 4.5 or above.

#### Animate Plugin
This plugin requires vtkFFMPEGIO module to be built in VTK.

#### Python Plugin
This plugin requires PythonQt 2.1 or above to be built AND installed.