/*=========================================================================
  Program: MILX MixView
  Module:
  Author: Jurgen Fripp
  Modified by: Shekhar Chandra
  Language: C++
  Created:

  Copyright: (c) 2009-11 CSIRO, Australia.

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.

  BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/

#ifndef __itkStatisticalShapeModelBase_h
#define __itkStatisticalShapeModelBase_h

// itk headers
#include <itkProcessObject.h>
#include <itkObjectFactory.h>
#include <itkMacro.h>
#include <itkPoint.h>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

// vtk headers
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataCollection.h>
#include <vtkFloatArray.h>
#include <vtkMatrix4x4.h>

// milx headers
#include "milxWin32Header.h"

namespace itk
{

/** 	\class StatisticalShapeModelBase
	\brief Base class for Statistical shape models. See class StatisticalShapeModel or RobustStatisticalShapeModel for more details.
 */

template <class TProfileSamplingPrecisionType=double>
class MILX_EXPORT StatisticalShapeModelBase : public ProcessObject
{
public:
  typedef StatisticalShapeModelBase     Self;
  typedef ProcessObject                 Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(StatisticalShapeModelBase, ProcessObject);

  typedef Point<double, 3> PointType;

  /**
   * Set Pose Type
   *  - 1, Rigid
   *  - 2, Similarity (Default)
   *  - 3, Affine
   *  - else none
   */
  itkSetClampMacro(PoseType, unsigned int, 1, 4)
  itkGetMacro( PoseType, unsigned int );

  /**
   * Load Model from file
   *  - Destroys current model
   */
  virtual bool LoadModel(const char *) throw(itk::ExceptionObject) { return false; }

  /**
   * Save Model as file
   *  - File format is as follows (float data)
   *  - 4 byte header (Number of Shapes (N), Dimensions (dims), Number of Points per shape (P), Number of Faces (F))
   *  - Data is written as (x_i,y_i,z_i) if dims == 3, where i = 1 to P for each shape 1... N
   * Therefore size should be 4 + dims*P*N + F
   */
  virtual bool SaveModel(const char *) { return false; }

  /**
   * Add shape to Model
   *  - Fileformat must be as outlined above
   *  - Warning vtkPoints * are currently Deleted in this function
   */
  bool AddShape(vtkPolyData *shape);
  bool SetShape(vtkPolyData *shape, unsigned int i);
  vtkPolyData* GetShape(unsigned int) throw(itk::ExceptionObject);
  /**
   * Get Number of Shapes in model
   */
  inline int GetNumberOfShapes()
  {
    return m_Landmarks.size();
  }
  /**
   * Get Number of Shapes in model
   */
  inline int GetNumberOfPoints()
  {
    return this->GetShape(0)->GetPoints()->GetNumberOfPoints();
  }

  /**
   * Remove shape from Model
   */
  virtual void RemoveShape(unsigned int);
  void RemoveAllShapes();
  /**
   *
   */
  void RemoveProcrustesAlignedPoints();

  /**
   * Get Shape After Procrustes Alignment
   */
  vtkPolyData * GetProcrustesAlignedSurface(int index);
  vtkPoints * GetAlignedPoints(unsigned int i)
    {
    this->Update();
    if (i >= (unsigned int)(m_Landmarks.size()))
      i = 0;
    //return m_ProcrustesAlign->GetOutput(i)->GetPoints();
    return m_ProcrustesAlignedPoints[i]->GetPoints();
    }
  bool AddAlignedShape(vtkPolyData *aligned);
  /**
   * Get Shape After Procrustes Alignment as is and quickly. This member will not trigger an Update of the filter.
   * Use GetProcrustesAlignedSurface() if worried about whether shape is up-to-date
   */
  vtkPolyData * GetAlignedShape(int index)
  {
    vtkPolyData * shape = vtkPolyData::New();

    shape->DeepCopy(GetMeanShape());
    shape->SetPoints(GetAlignedPoints(index));

    return shape;
  }

  /**
   * Get Centroid
   */
  PointType GetCentroid(vtkPolyData * surface);

  /**
   * Get Centroid Size
   */
  double GetCentroidSize(vtkPolyData * surface);
  double GetCentroidSize(vtkPolyData * surface, PointType centroid);

  /**
   * Set Mode precision of PCA
   */
  itkSetClampMacro(Precision, float, 0, 1.0);
  itkGetMacro(Precision, float);
  virtual int GetNumberOfModes() { return 0; }
  virtual float GetPrecisionRequiredFrom(int modes) { return 0.9; }

  /**
   * Get the surface described by
   *  - translation t
   *  - scale and rotation (theta, phi, psi) where [0: 2pi] [0 : pi] [0 : 2pi]
   *  - b shape parameters
   *  - centroid (if true matrix doesn't have centroid info included in it)
   */
  virtual vtkPolyData * GetSurface(double s, PointType t, PointType orientation, vtkFloatArray *b = NULL);
  virtual vtkPolyData * GetSurface(double s, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b = NULL);
  virtual vtkPolyData * GetSurface(vnl_matrix<double> pose, vtkFloatArray *b = NULL, bool centroid = true);
  virtual vtkPolyData * GetSurface(vtkMatrix4x4 * pose, vtkFloatArray *b = NULL, bool centroid = true);

  /**
   * Get parameters that describe surface
   * If useParams true then fit uses t, s, orientation shape to compare not the mean shape
   */
  virtual void GetSurfaceSimilarityParameters(vtkPolyData * surface, double &s, PointType &t, PointType &orientation, vtkFloatArray *b, float bConstrain=3, unsigned int maxModesCount=0);
  virtual void GetSurfaceSimilarityParameters(vtkPolyData * surface, double &s, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b, float bConstrain=3, unsigned int maxModesCount=0);
  virtual void GetSurfaceSimilarityParameters(vtkPolyData * surface, vtkFloatArray *weights, vtkMatrix4x4 *matrix, vtkFloatArray *b) {}

  virtual void GetSurfacePose(vtkPolyData * surface, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b = NULL);
  virtual void GetSurfacePose2(vtkPolyData * surface, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b = NULL);
  virtual void GetSurfaceShapeParams(vtkPolyData * surface, double scale, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b);

  /**
   * Check whether the difference (p2p) is less than threshold
   */
  bool CheckConvergence(vtkPolyData * surface1, vtkPolyData * surface2, float minimumAverageSquareDistance = 0.001);

  /**
   * Calculate the square distance between two surfaces
   *  - Assumes polyData has same number of points
   */
  double PolyDataDistance(vtkPolyData * poly1, vtkPolyData * poly2, vtkFloatArray* errorArray = NULL);
  double PolyDataSquareDistance(vtkPolyData * poly1, vtkPolyData * poly2, vtkFloatArray* errorArray = NULL);
  /**
   * Calculates the Mahalanobis distance between two surfaces within the shape model
   *  - Assumes polyData has same number of points
   *  - Assumes the surfaces have been pre co-registered.
   */
  virtual double MahalanobisDistanceBetween(vtkPolyData * poly1, vtkPolyData * poly2, const bool includeScaling = false) { return 0.0; }

  /**
   * Generate Shape using b
   */
  virtual vtkPolyData * GetParameterisedShape(vtkFloatArray *b) { return NULL; }
  /**
   * Generate Parameters that best describe shape
   */
  virtual vtkFloatArray * GetShapeParameters(vtkPolyData *shape, int bsize = -1) { return NULL; }

  /**
   * Get Mean Shape
   */
  vtkPolyData * GetMeanShape();
  double GetMeanShapeSize()
    {
    this->Update();
    return m_MeanScale;
    }

  virtual void Update();

  /**
   * Clip all parts of the existing surface based on its scalars. Parts to be clipped are assumed marked as 0.0 and kept parts as 1.0
   * One gets the id mapping for the missing data, so that correspondence can be maintained.
   * Usage:
   \code
   vtkSmartPointer<vtkPolyData> clippedSurfaceX = vtkSmartPointer<vtkPolyData>::New();
   vtkSmartPointer<vtkIdTypeArray> idMap = vtkSmartPointer<vtkIdTypeArray>::New();
   ClipSurfaceBasedOnScalars(surfaceX, clippedSurfaceX, idMap);
   \endcode
   */
  static void ClipSurfaceBasedOnScalars(vtkPolyData *surface, vtkPolyData *clippedSurface, vtkIdTypeArray *ids, const float thresholdValue = 1.0);
  /**
   * Restore correspondence of clipped mesh based on point id map provided. The empty parts are filled with the point corresponding to point id 0.
   *
   */
  static void RestoreClippedSurfaceCorrespondence(vtkPolyData *clippedSurface, vtkIdTypeArray *ids, vtkPolyData *clippedCorrespondSurface);

  //Robust and Weighted virtual members
  /**
   * Set/Get whether to use a shape constraints based on the Weighted SSM
   */
  itkSetMacro( WeightMode, bool );
  itkGetMacro( WeightMode, bool );
  itkBooleanMacro(WeightMode);
  /**
   * Set/Get whether to use a shape constraints based on the Robust SSM
   */
  itkSetMacro( RobustMode, bool );
  itkGetMacro( RobustMode, bool );
  itkBooleanMacro(RobustMode);
  /**
   * Set/Get whether to use a shape constraints based on the Standard SSM. Default
   */
  itkSetMacro( StandardMode, bool );
  itkGetMacro( StandardMode, bool );
  itkBooleanMacro(StandardMode);

  virtual bool GetMissingMode() { return false; }

protected:
  StatisticalShapeModelBase();
  virtual ~StatisticalShapeModelBase();

  virtual void GenerateData() {}

  /**
   * Get parameters that describe surface
   * If useParams true then fit uses t, s, orientation shape to compare not the mean shape
   */
  void GetSurfaceParametersToMatchShape1(vtkPolyData * surface, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b, float bConstrain=3, unsigned int maxModesCount=0);

  /**
   * Update Procrustes Alignment and PCA
   */
  virtual void UpdateProcrustesAlignment() throw(itk::ExceptionObject);
  virtual void UpdatePCA() throw(itk::ExceptionObject) {}

  /**
   * Set/Get the Valid Data Bit for the class
   *  - Used by Update
   */
  itkSetMacro( Valid, bool );
  itkGetMacro( Valid, bool );

  /**
   * Generate Similarity matrix from orientation, translation and scale (around centroid)
   */
  virtual vtkMatrix4x4 *GetSimilarityMatrix(double s, double tx, double ty, double tz, double theta, double phi, double psi, double centroidX = 0, double centroidY = 0, double centroidZ = 0);
  virtual vtkMatrix4x4 *GetSimilarityMatrix(double s, PointType t, PointType orientation, PointType centroid);
  virtual void GetSimilarityMatrix(vnl_matrix<double> &matrix, double s, PointType t, PointType orientation, PointType centroid);

  /**
   * Get Euler angles from Rotation matrix
   */
  void GetSurfacePoseFromTransformMatrix(vtkMatrix4x4 * matrix, PointType &orientation);
  void GetSurfacePoseFromTransformMatrix(vnl_matrix<double> matrix, PointType &orientation);

  //
  bool                          m_Valid;
  // Store as polydata as needed as filter input in this form
  std::vector<vtkPolyData *>    m_Landmarks;
  // The below polyData has no surface information
  std::vector<vtkPolyData *>    m_ProcrustesAlignedPoints;
  // Mean shape
  vtkPolyData *                 m_MeanShape;

  double       m_MeanScale;
  unsigned int m_PoseType;
  int          m_Modes;
  float        m_Precision;

  bool m_WeightMode; //!< Weighted shape model mode
  bool m_RobustMode; //!< Robust shape model mode
  bool m_StandardMode; //!< Standard shape model mode

private:
  StatisticalShapeModelBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStatisticalShapeModelBase.txx"
#endif

#endif
