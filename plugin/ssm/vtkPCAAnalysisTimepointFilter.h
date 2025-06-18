/*=========================================================================
  Program: MILX MSK
  Module: milx-vtkHybrid
  Author: Shekhar Chandra
  Modified by:
  Language: C++
  Created: Tue 19/04/2011

  Copyright: (c) 2009-12 CSIRO, Australia.

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.

  BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/

#ifndef __vtkPCAAnalysisTimepointFilter_h
#define __vtkPCAAnalysisTimepointFilter_h

#include <vtkFloatArray.h>
#include <vtkPCAAnalysisFilter.h>

#include "milxWin32Header.h"

class vtkPointSet;

class MILX_EXPORT vtkPCAAnalysisTimepointFilter : public vtkPCAAnalysisFilter
{
 public:
  vtkTypeMacro(vtkPCAAnalysisTimepointFilter,vtkPCAAnalysisFilter);

  // Description:
  // Prints information about the state of the filter.
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Creates with similarity transform.
  static vtkPCAAnalysisTimepointFilter *New();

  // Description:
  // Get the vector of eigenvalues sorted in descending order
  vtkSetObjectMacro(Evals, vtkFloatArray);
  vtkGetObjectMacro(Evals, vtkFloatArray);

  void SetEvectors(double** vecs)
  {
    evecMat2 = vecs;
  }
  double** GetEvectors()
  {
    return evecMat2;
  }

  void SetMeanShape(double* shape)
  {
    meanshape = shape;
  }
  double* GetMeanShape()
  {
    return meanshape;
  }

  // Description:
  // Timepoint mode is to replace the mean shape with the first shape. Useful to make a model out of timepoint shapes
  vtkSetMacro(TimepointMode, bool);
  vtkGetMacro(TimepointMode, bool);
  vtkBooleanMacro(TimepointMode, bool);

  //!Get/Set whether eigenvalues and eigenvectors are already known or pre-loaded, if so, don't bother doing EM algorithm
  //!Use stored result.
  vtkSetMacro(PreLoaded, bool);
  vtkGetMacro(PreLoaded, bool);
  vtkBooleanMacro(PreLoaded, bool);

  // Description:
  // Fills the shape with:
  //
  // mean + b[0] * sqrt(eigenvalue[0]) * eigenvector[0]
  //      + b[1] * sqrt(eigenvalue[1]) * eigenvector[1]
  // ...
  //      + b[sizeb-1] * sqrt(eigenvalue[bsize-1]) * eigenvector[bsize-1]
  //
  // here b are the parameters expressed in standard deviations
  // bsize is the number of parameters in the b vector
  // This function assumes that shape is allready allocated
  // with the right size, it just moves the points.
  void GetParameterisedShape(vtkFloatArray *b, vtkPointSet* shape);

  // Description:
  // Return the bsize parameters b that best model the given shape
  // (in standard deviations).
  // That is that the given shape will be approximated by:
  //
  // shape ~ mean + b[0] * sqrt(eigenvalue[0]) * eigenvector[0]
  //              + b[1] * sqrt(eigenvalue[1]) * eigenvector[1]
  //         ...
  //              + b[bsize-1] * sqrt(eigenvalue[bsize-1]) * eigenvector[bsize-1]
  void GetShapeParameters(vtkPointSet *shape, vtkFloatArray *b, int bsize);

  // Description:
  // Retrieve how many modes are necessary to model the given proportion of the variation.
  // proportion should be between 0 and 1
  int GetModesRequiredFor(double proportion);

protected:
  vtkPCAAnalysisTimepointFilter();
  virtual ~vtkPCAAnalysisTimepointFilter();

  // Description:
  // Usual data generation method.
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // Eigenvalues
  vtkFloatArray *Evals;

  // Matrix where each column is an eigenvector
  double **evecMat2;

  // The mean shape in a vector
  double *meanshape;

  bool TimepointMode; //!< Timepoint enabled?
  bool PreLoaded; //!< Pre loaded all data?

private:
  vtkPCAAnalysisTimepointFilter(const vtkPCAAnalysisTimepointFilter&);  // Not implemented.
  void operator=(const vtkPCAAnalysisTimepointFilter&);  // Not implemented.

};
#endif
