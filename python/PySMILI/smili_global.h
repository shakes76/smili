/*
  This file is part of SMILI.
*/
#ifndef PYTHON_BINDINGS 
#define PYTHON_BINDINGS

#pragma once

// Make "signals:", "slots:" visible as access specifiers
#define QT_ANNOTATE_ACCESS_SPECIFIER(a) __attribute__((annotate(#a)))

//ITK
//#include <itkImage.h>

//VTK
//#include <vtkSmartPointer.h>
//#include <vtkPointSet.h> //base of vtkPolyData
#include <vtkPolyData.h>
//#include <vtkImageData.h>
#if VTK_MAJOR_VERSION <= 8
  #include <QVTKWidget.h>
#else
  #include <QVTKOpenGLStereoWidget.h>
#endif

#if VTK_MAJOR_VERSION > 8
  typedef QVTKOpenGLStereoWidget QVTKWidget;
#endif

//SMILI
//#include "milxMath.h"
// #include "milxModel.h"
// #include "milxFile.h"
//milxQt
#include "milxQtWindow.h"
#include "milxQtUnifiedWindow.h"
#include "milxQtRenderWindow.h"
#include "milxQtModel.h"
#include "milxQtImage.h"
#include "milxQtFile.h"
#include "milxQtPlot.h"
// #include "milxQtMain.h"

#endif // Define PYTHON_BINDINGS
