/*
  This file is part of SMILI.
*/
#ifndef PYTHON_BINDINGS
#define PYTHON_BINDINGS

#pragma once

// Make "signals:", "slots:" visible as access specifiers
#define QT_ANNOTATE_ACCESS_SPECIFIER(a) __attribute__((annotate(#a)))

//ITK
//#include <itkArray.h>
//#include <itkImage.h>

//VTK
// #include <vtkSmartPointer.h>
// #include <vtkMatrix4x4.h>
// #include <vtkDataObject.h>
// #include <vtkObjectBase.h>s
// #include <vtkObject.h>
#include <vtkPoints.h>
#include <vtkDataObject.h> //base of vtkDataSet
#include <vtkPointSet.h> //base of vtkPolyData
#include <vtkPolyData.h>
// #include <vtkDataSet.h> //base of vtkImageData
#include <vtkImageData.h>
#if VTK_MAJOR_VERSION <= 8
  #include <QVTKWidget.h>
#else
  #include <QVTKOpenGLNativeWidget.h>
#endif

//SMILI
//#include "milxMath.h"
// #include "milxModel.h"
// #include "milxFile.h"
//milxQt
#include "milxQtWindow.h"
//#include "milxQtUnifiedWindow.h"
#include "milxQtRenderWindow.h"
#include "milxQtModel.h"
#include "milxQtImage.h"
#include "milxQtFile.h"
#include "milxQtPlot.h"
#include "milxQtMain.h"

#endif // Define PYTHON_BINDINGS
