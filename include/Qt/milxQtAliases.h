/*=========================================================================
  The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO)
  ABN 41 687 119 230.
  All rights reserved.

  Licensed under the CSIRO BSD 3-Clause License
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License in the file LICENSE.md or at

  https://stash.csiro.au/projects/SMILI/repos/smili/browse/license.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
=========================================================================*/
/**
* \file milxQtAliases.h
*
* \brief This file defines all the defines, aliases and frequently used functions and variables.
* \author Shekhar S. Chandra, 2013.
*
* This file is part of MILX Qt Library.
*
*/
#ifndef MILXQTALIASES_H_INCLUDED
#define MILXQTALIASES_H_INCLUDED
/**
* \mainpage MILX Qt API
*
* \section intro_sec Introduction
*
* The MILX Qt library essentially binds together three important libraries.
* Qt GUI Library for user interfaces. The VTK Visualisation library for OpenGL data plot rendering, and the
* Insight Toolkit library for image Registration and Segmentation.
*
* The library consists of four main classes intended for direct use. These are milxQtRenderWindow, milxQtImage, milxQtModel, milxQtFile and milxQtMain.
* milxQtMain is used as the main GUI class and is able to load plugins. Plugins are created using the milxQtPluginInterface class.
* The main GUI application is called sMILX and uses the SMILI library for all image and model processing.
* See the plugin guide for how to create plugins for sMILX.
*/
#include <limits>
#include <stdexcept> ///Include stdexcept to handle conversion exceptions
#include <sstream>
#include <iostream>
#include <vnl/vnl_vector_fixed.h>
#include "milxGlobal.h"

using namespace std;
using milx::coordinate; //allow to enable cleaner looking code
using milx::coordinateType; //allow to enable cleaner looking code

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif

///Define Windows DLL importing
#if defined(WIN32)
#ifdef MILXQT_DLL
/*
Used for DLL generation purposes (Windows Specific) Import/Export.
Templates classes cannot be imported hence its own variable.
The export command is
*/
#if defined(MILXQT_MAKEDLL)     // create a MILXQT DLL library
#define MILXQT_EXPORT  __declspec(dllexport)
#define MILXQT_TEMPLATEDLL  __declspec(dllexport)
#else                        // use a MILXQT DLL library
#define MILXQT_EXPORT  __declspec(dllimport)
#define MILXQT_TEMPLATEDLL
#endif
#endif // MILXQT_DLL
#endif // WIN32

#ifndef MILXQT_EXPORT
#define MILXQT_EXPORT
#define MILXQT_TEMPLATEDLL
#endif

///Define Windows Plugin DLL importing
#if defined(WIN32)
#ifdef MILXQT_PLUGIN_DLL
/*
Used for DLL generation purposes (Windows Specific) Import/Export.
Templates classes cannot be imported hence its own variable.
The export command is
*/
#if defined(MILXQT_PLUGIN_MAKEDLL)     // create a MILXQT DLL library
#define MILXQT_PLUGIN_EXPORT  __declspec(dllexport)
#define MILXQT_PLUGIN_TEMPLATEDLL  __declspec(dllexport)
#else                        // use a MILXQT DLL library
#define MILXQT_PLUGIN_EXPORT  __declspec(dllimport)
#define MILXQT_PLUGIN_TEMPLATEDLL
#endif
#endif // MILXQT_PLUGIN_DLL
#endif // WIN32

#ifndef MILXQT_PLUGIN_EXPORT
#define MILXQT_PLUGIN_EXPORT
#define MILXQT_PLUGIN_TEMPLATEDLL
#endif

static const float milxQtVersion = static_cast<float>(1.04);
static const int minWindowSize = 256;
static const int maxAASamples = 2; //Anti-Aliasing

#ifndef DEF_EXTS
#define DEF_EXTS
/**
    These common constant strings are used in file open/save dialogs
*/
//Opens
static string openImageExts = "Images (*.png *.jpeg *.jpg *.bmp *.tiff *.tif *.pbm *.pgm *.ppm)";
static string openSeriesImageExts = "Image Series (*.ima *.dcm *.dicom *.gz)";
static string openMedImageExts = "Medical Images (*.nii *.gz *.ima *.dcm *.dicom *.mha *.mhd *.img *.hdr *.gipl *.spr)";
static string openOtherImageExts = "Other Images (*.vti)"; ///\todo add raw support, not working atm
static string openOtherExts = "Text Delimited Files (*.csv *.dat *.txt)";
static string openModelExts = "Model or Polygonal Files (*.vtp *.vtk *.ply *.obj *.stl)";
static string extensionsOpen = "*.png *.jpeg *.jpg *.bmp *.tiff *.tif *.pbm *.pgm *.ppm *.nii *.gz *.ima *.dcm *.dicom *.mha *.mhd *.img *.hdr *.gipl *.spr *.vti *.vtp *.vtk *.ply *.obj *.stl *.csv *.dat *.txt";
static string allFileExts = "All Files (*.*)";
static string openSupportedExts = "Images and Model Files (*.png *.jpeg *.jpg *.bmp *.tiff *.tif *.pbm *.pgm *.ppm *.nii *.gz *.ima *.dcm *.dicom *.mhd *.img *.hdr *.vti *.mrc *.rec *.vtp *.vtk *.ply *.obj *.stl)";
static string openExts = openSupportedExts + ";;" + openModelExts + ";;" + openMedImageExts + ";;" + openImageExts + ";;" + openOtherImageExts + ";;" + openOtherExts + ";;" + allFileExts;
//Saves
static string saveImageExts = "Images (*.png *.jpeg *.jpg *.bmp *.tiff *.tif)";
static string saveMedImageExts = "Medical Images (*.nii *.gz *.ima *.dcm *.dicom *.mha *.mhd *.img *.hdr)";
static string saveOtherExts = "Text Delimited Files (*.csv *.dat *.txt)";
static string saveOtherImageExts = "Other Images (*.vti)";
static string saveModelExts = "Model or Polygonal File (*.vtp *.vtk *.ply *.obj *.stl)";
static string extensionsSave = "*.png *.jpeg *.jpg *.bmp *.tiff *.tif *.nii *.gz *.ima *.dcm *.dicom *.mhd *.img *.hdr *.vti *.vtp *.vtk *.ply *.obj *.stl";
static string saveExtsForImages = saveMedImageExts + ";;" + saveImageExts + ";;" + saveOtherImageExts;
static string saveExtsForScreens = saveImageExts;
static string saveExtsForModels = saveModelExts;
#endif //DEF_EXTS

#endif // MILXQTALIASES_H_INCLUDED
