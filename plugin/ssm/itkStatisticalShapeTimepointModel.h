/*=========================================================================
  Program: MILX MixView
  Module: milxSSM3D
  Author: Shekhar Chandra
  Language: C++
  Created:

  Copyright: (c) 2009-12 CSIRO, Australia.

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.

  BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/

#ifndef __itkStatisticalShapeTimepointModel_h
#define __itkStatisticalShapeTimepointModel_h

// milx headers
#include "milxWin32Header.h"
#include "itkStatisticalShapeModel.h"
#include "vtkPCAAnalysisTimepointFilter.h"

namespace itk
{

template <class TProfileSamplingPrecisionType=double>
class MILX_EXPORT StatisticalShapeTimepointModel : public StatisticalShapeModel<TProfileSamplingPrecisionType>
{
public:
  typedef StatisticalShapeTimepointModel                                  Self;
  typedef StatisticalShapeModel<TProfileSamplingPrecisionType>            Superclass;
  typedef SmartPointer<Self>                                              Pointer;
  typedef SmartPointer<const Self>                                        ConstPointer;

  typedef Point<double, 3> PointType;

  itkTypeMacro(StatisticalShapeTimepointModel, StatisticalShapeModel<TProfileSamplingPrecisionType>);

  itkNewMacro(Self);

  bool AddAlignedShape(vtkPolyData *aligned);

  /**
   * Load Model from file fully including eigenvalues and eigenvectors
   *  - Destroys current model
   */
  virtual bool LoadModel(const char *) throw(itk::ExceptionObject);

  /**
   * Save Model as file fully including eigenvalues and eigenvectors
   *  - File format is as follows (float data)
   *  - 4 byte header (Number of Shapes (N), Dimensions (dims), Number of Points per shape (P), Number of Faces (F))
   *  - Data is written as (x_i,y_i,z_i) if dims == 3, where i = 1 to P for each shape 1... N
   * Therefore size should be 4 + dims*P*N + F
   */
  virtual bool SaveModel(const char *);

  /**
   * Remove shape from Model
   */
  virtual void RemoveShape(unsigned int);

  /**
   * Generate Shape using b
   */
  virtual vtkPolyData * GetParameterisedShape(vtkFloatArray *b);

  /**
   * Generate Parameters that best describe shape
   */
  vtkFloatArray * GetShapeParameters(vtkPolyData *shape, int bsize = -1);

  //Set/Get Timepoint mode
  itkSetMacro(TimepointMode, bool);
  itkGetMacro(TimepointMode, bool);
  itkBooleanMacro(TimepointMode);

  virtual void Update();

  virtual int GetNumberOfModes()
  {
    this->Update();
    return this->GetPCA()->GetModesRequiredFor(this->m_Precision);
  }
  virtual float GetPrecisionRequiredFrom(int modes);

  virtual vtkPolyData * GetSurface(double s, PointType t, PointType orientation, vtkFloatArray *b = NULL);
  virtual vtkPolyData * GetSurface(double s, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b = NULL);
  virtual vtkPolyData * GetSurface(vnl_matrix<double> pose, vtkFloatArray *b = NULL, bool centroid = true);
  virtual vtkPolyData * GetSurface(vtkMatrix4x4 * pose, vtkFloatArray *b = NULL, bool centroid = true);
  virtual void GetSurfacePose(vtkPolyData * surface, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b = NULL);
  virtual void GetSurfacePose2(vtkPolyData * surface, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b = NULL);
  virtual void GetSurfaceShapeParams(vtkPolyData * surface, double scale, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b);

  void Reset();

  /**
    * Allow access to internal structures
    *  - Basically only needed for debugging and visualization (ie should be hidden or removed later)
    */
  vtkPCAAnalysisTimepointFilter* GetPCA();

protected:
  StatisticalShapeTimepointModel();
  virtual ~StatisticalShapeTimepointModel();

  virtual void GenerateData();

  /**
   * Update Procrustes Alignment and PCA
   */
  virtual void UpdateProcrustesAlignment() throw(itk::ExceptionObject);
  virtual void UpdatePCA() throw(itk::ExceptionObject);

  /**
   * Generate Similarity matrix from orientation, translation and scale (around centroid)
   */
  virtual vtkMatrix4x4 *GetSimilarityMatrix(double s, double tx, double ty, double tz, double theta, double phi, double psi, double centroidX = 0, double centroidY = 0, double centroidZ = 0);
  virtual vtkMatrix4x4 *GetSimilarityMatrix(double s, PointType t, PointType orientation, PointType centroid);
  virtual void GetSimilarityMatrix(vnl_matrix<double> &matrix, double s, PointType t, PointType orientation, PointType centroid);

  bool m_TimepointMode; //!< Timepoint enabled?
  bool m_Loaded; //!< Loaded from file?
  bool m_PreAligned; //!< pre aligned shapes?

private:
  StatisticalShapeTimepointModel(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStatisticalShapeTimepointModel.txx"
#endif

#endif
