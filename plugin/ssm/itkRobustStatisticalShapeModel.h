/*=========================================================================
  Program: MILX MixView
  Module:
  Author:
  Modified by:
  Language: C++
  Created:

  Copyright: (c) 2009 CSIRO, Australia.

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.

  BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/

#ifndef __itkRobustStatisticalShapeModel_h
#define __itkRobustStatisticalShapeModel_h

// vtk headers
#include <vtkWeightedPCAAnalysisFilter.h>

#include "itkStatisticalShapeModelBase.h"

namespace itk
{

/** 	\class RobustStatisticalShapeModel
	\brief Construct, save and open robust Statistical Shape Models (SSMs).

	The shapes passed (via AddShape()) are aligned using Procrustes Alignment and eigen decomposed via a Weighted PCA.
  The weights can be either 0.0 for missing points or 1.0 for known points.
	Shapes can then be parameterised using the principal components of the RSSM.

	An example application includes computing the model of the hip bone from missing data because of the small field of view.

  The weights must be set via a float array passed also within AddShape. The model can be loaded from file using
  LoadRobustModel and saved via SaveRobustModel. These members write the eigen-decomposition as well as the meshes used for training.
  Members LoadModel and SaveModel don't do this.

  Binary Weighted Procrustes is computed. All regions with scalars below 1e-3 are ignored in alignment.

  Special functions are provided for clipping parts of surfaces (such as being outside an image) depending on their scalar values. These assume the scalars of the mesh are 0.0 for clipping and 1.0 for keeping points.
  These are ClipSurfaceOutsideImage() and ClippedSurfaceWithCorrespondence().

  Notes: Shapes are not deleted, unless RemoveAllShapes() is invoked by user. Alignment can be disabled by PreAlignedOn(), which will assume shapes are already aligned.
  Also, weights are not automatically set as mesh scalars, so when projecting shapes (even training shapes), scalars need to be set manually by you.

	Declaration Example:
	\code
	typedef itk::RobustStatisticalShapeModel<double> ShapeModelType;

	ShapeModelType::Pointer m_SSM; //!< The Statistical Shape Model
	m_SSM = ShapeModelType::New();
	\endcode

	Usage Example 1:
	\code
	vtkPolyDataCollection* meshes = loadMeshes(); //Load some of your meshes here, implement yourself

	size_t n = meshes->GetNumberOfItems();
	meshes->InitTraversal();
	for(size_t j = 0; j < n; j ++)
  {
    //Determine weights....
    weights = ...

		m_SSM->AddShape(meshes->GetNextItem(), weights); //!< Set model to shape filter
  }

	m_SSM->Update();
	//Dont destroy meshes collection, SSM class handles that
	\endcode

	Usage Example 2 Load Robust Model:
	\code
	m_SSM->LoadRobustModel(filename.toStdString().c_str());
	m_SSM->Update();
	\endcode

  Usage Example 3 Load Normal SSM Model File:
	\code
	m_SSM->LoadModel(filename.toStdString().c_str());
	m_SSM->Update();
	\endcode

  Further Notes:
  To project mesh and recover similarity parameters, use the matrix method and not the traditional variable method (scale, tx, ty, etc.):
  \code
  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
  ssm->GetSurfaceSimilarityParameters(surface, weights, matrix, b);
  vtkSmartPointer<vtkPolyData> constrainedSurface = ssm->GetSurface(matrix, b, false);
  \endcode
 */

template <class TProfileSamplingPrecisionType=double>
class MILX_EXPORT RobustStatisticalShapeModel : public StatisticalShapeModelBase<TProfileSamplingPrecisionType>
{
public:
  typedef RobustStatisticalShapeModel                                   Self;
  typedef StatisticalShapeModelBase<TProfileSamplingPrecisionType>      Superclass;
  typedef SmartPointer<Self>                                            Pointer;
  typedef SmartPointer<const Self>                                      ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(RobustStatisticalShapeModel, StatisticalShapeModelBase);

  typedef Point<double, 3> PointType;

  /**
   * Load Model from file
   *  - Destroys current model
   */
  virtual bool LoadModel(const char *) throw(itk::ExceptionObject);
  virtual bool LoadCompactModel(const char *) throw(itk::ExceptionObject);

  /**
   * Save Model as file
   *  - File format is as follows (float data)
   *  - 4 byte header (Number of Shapes (N), Dimensions (dims), Number of Points per shape (P), Number of Faces (F))
   *  - Data is written as (x_i,y_i,z_i) if dims == 3, where i = 1 to P for each shape 1... N
   * Therefore size should be 4 + dims*P*N + F
   */
  virtual bool SaveModel(const char *);
  virtual bool SaveCompactModel(const char *);

  /**
   * Add shape to Model
   *  - Fileformat must be as outlined above
   *  - Warning vtkPoints * are currently Deleted in this function
   */
  virtual bool AddShape(vtkPolyData *shape, vtkFloatArray *weights);
  virtual bool SetShape(vtkPolyData *shape, vtkFloatArray *weights, unsigned int i);
  ///Return the clipped weighted surface, the lowest scalar value will be used as threshold below
  vtkPolyData* GetWeightedShape(unsigned int) throw(itk::ExceptionObject);
  /**
   * Remove shape from Model
   */
  virtual void RemoveShape(unsigned int);

  /**
   * Remove shape from Model
   */
  virtual void RemoveAllWeights();

  PointType GetWeightedCentroid(vtkPolyData * surface, vtkFloatArray *weights);
  double GetWeightedCentroidSize(vtkPolyData * surface, vtkFloatArray *weights, PointType centroid);

  virtual int GetNumberOfModes()
  {
    this->Update();
    return this->GetPCA()->GetModesRequiredFor(this->m_Precision);
  }
  virtual float GetPrecisionRequiredFrom(int modes);

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
   * Get parameters that describe surface that is in real space
   * If useParams true then fit uses t, s, orientation shape to compare not the mean shape
   * \code
    vtkSmartPointer<vtkFloatArray> b = vtkSmartPointer<vtkFloatArray>::New();
      b->SetNumberOfTuples(ssm->GetNumberOfModes());
      b->SetNumberOfComponents(1);
      b->FillComponent(0, 0.0);
    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    ssm->GetSurfaceSimilarityParameters(mesh, weights, matrix, b);
    vtkSmartPointer<vtkPolyData> constrainedSurface = ssm->GetSurface(matrix, b, false);
   * \endcode
   */
  virtual void GetSurfaceSimilarityParameters(vtkPolyData * surface, vtkFloatArray *weights, vtkMatrix4x4 *matrix, vtkFloatArray *b);

  virtual void GetSurfacePose(vtkPolyData * surface, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b = NULL);
  virtual void GetWeightedSurfacePose(vtkPolyData * surface, vtkFloatArray *weights, vtkFloatArray *b = NULL);
  virtual void GetMissingSurfacePose2(vtkPolyData * surface, vtkFloatArray *weights, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b = NULL);
  virtual void GetMissingSurfacePose(vtkPolyData * surface, vtkFloatArray *weights, vtkFloatArray *b = NULL);
  virtual void GetSurfacePose2(vtkPolyData * surface, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b = NULL);
  virtual void GetSurfaceShapeParams(vtkPolyData * surface, vtkFloatArray *weights, double scale, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b);
  virtual void GetWeightedSurfaceShapeParams(vtkPolyData * surface, vtkFloatArray *weights, vtkFloatArray *b);
  virtual void GetMissingSurfaceShapeParams2(vtkPolyData * surface, vtkFloatArray *weights, double scale, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b);
  virtual void GetMissingSurfaceShapeParams(vtkPolyData * surface, vtkFloatArray *weights, vtkFloatArray *b);

  /**
   * Calculates the Mahalanobis distance between two surfaces within the shape model
   *  - Assumes polyData has same number of points
   *  - Assumes the surfaces have been pre co-registered.
   */
  virtual double MahalanobisDistanceBetween(vtkPolyData * poly1, vtkFloatArray *weights1, vtkPolyData * poly2, vtkFloatArray *weights2);

  /**
   * Generate Shape using b
   */
  virtual vtkPolyData * GetParameterisedShape(vtkFloatArray *b);
  virtual vtkPolyData * GetParameterisedShape(vtkFloatArray *b, vtkFloatArray *weights);
  /**
   * Generate Parameters that best describe shape
   */
  virtual vtkFloatArray * GetShapeParameters(vtkPolyData *shape, vtkFloatArray *weights, int bsize = -1);
  virtual vtkFloatArray * GetShapeParametersInPlace(vtkPolyData *shape, vtkFloatArray *weights, int bsize = -1);

  /**
    * Allow access to internal structures
    *  - Basically only needed for debugging and visualization (ie should be hidden or removed later)
    */
  inline vtkWeightedPCAAnalysisFilter* GetPCA() { return m_PCA; }

  //Robust methods
  /** \brief Get the weights, only a pointer pass.
   *
   * It is assumed that the number of points are the same for each shape.
   * \return vtkFloatArray*
   *
   */
  vtkFloatArray* GetWeights(unsigned int) throw(itk::ExceptionObject);

  /**
   * Tolerance for checking convergence for the EM algorithm, default is 1e-9
   */
  itkSetMacro(Tolerance, float)
  itkGetMacro(Tolerance, float)

  /**
   * Max iterations of the EM algorithm, default is 3
   */
  itkSetMacro(MaxIterations, int)
  itkGetMacro(MaxIterations, int)

  /**
   * Set/Get the whether the shapes are already aligned.
   */
  itkSetMacro( PreAligned, bool );
  itkGetMacro( PreAligned, bool );
  itkBooleanMacro( PreAligned );

  /**
   * Set/Get the whether the shapes are scale invariant. If true, all alignment and reconstruction is done
   * using unit scale, i.e. scale is always stripped. If set to false, you must ensure that your training surfaces do not have their scale stripped also.
   */
  itkSetMacro( ScaleInvariance, bool );
  itkGetMacro( ScaleInvariance, bool );
  itkBooleanMacro( ScaleInvariance );

  /**
   * Set/Get the loaded shape model solution should be ignored.
   * This can be used to override the loaded shape model solution by forcing
   * recomputation of the RPCA.
   */
  itkSetMacro( Loaded, bool );
  itkGetMacro( Loaded, bool );
  itkBooleanMacro( Loaded );

  /**
   * Set/Get the whether to use random guess as initialisation for the WPCA
   * Only applicable for weighted mode. Default false.
   * If false, a perturbed PCA solution is used to initialise the WPCA.
   */
  itkSetMacro( RandomInitialisation, bool );
  itkGetMacro( RandomInitialisation, bool );
  itkBooleanMacro( RandomInitialisation );

  /**
   * Set/Get the mode of the shape model solution.
   * This can be used to toggle between weighted and missing modes within
   * the RPCA. (Default missing mode.). You must ensure the weights are set accordingly for the relevant mode.
   */
  itkSetMacro( MissingMode, bool );
  itkGetMacro( MissingMode, bool );
  itkBooleanMacro( MissingMode );
  /**
   * Set/Get the mode of the shape model solution.
   * use this to fine tune the reconstructions given missing parts.
   */
  inline void SetMissingReconstuctionIterations(const size_t iterations)
  { m_PCA->SetReconstructionIterations(iterations);  }
  inline size_t GetMissingReconstuctionIterations()
  { return m_PCA->GetReconstructionIterations(); }

  /**
   * Set Mode precision of Robust PCA Iterative Reconstructions
   */
  itkSetClampMacro(RobustPCAPrecision, float, 0, 1.0)
  itkGetMacro(RobustPCAPrecision, float)

protected:
  RobustStatisticalShapeModel();
  virtual ~RobustStatisticalShapeModel();

  float m_Tolerance;
  int m_MaxIterations;

  virtual void GenerateData();

  /**
   * Get parameters that describe surface
   * If useParams true then fit uses t, s, orientation shape to compare not the mean shape
   */
  virtual void GetSurfaceParametersToMatchShape1(vtkPolyData * surface, vtkFloatArray *weights, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b);
  virtual void GetSurfaceParametersToMatchShape2(vtkPolyData * surface, vtkFloatArray *weights, vtkMatrix4x4 *matrix, vtkFloatArray *b);
  virtual void GetSurfaceParametersToMatchShape3(vtkPolyData * surface, vtkFloatArray *weights, vtkMatrix4x4 *matrix, vtkFloatArray *b);

  /**
   * Update Procrustes Alignment and PCA
   */
  virtual void UpdateProcrustesAlignment() throw(itk::ExceptionObject);
  virtual void UpdatePCA() throw(itk::ExceptionObject);

  //Not in use so disable
  virtual void GetSurfaceSimilarityParameters(vtkPolyData * surface, vtkFloatArray *weights, double &s, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b);

  //
  bool m_Valid;
  bool m_Loaded;
  bool m_PreAligned;
  bool m_ScaleInvariance;
  bool m_MissingMode;
  bool m_RandomInitialisation;

  // PCA
  vtkWeightedPCAAnalysisFilter* m_PCA;
  //Weights for the points in each shape
  std::vector<vtkFloatArray *>  m_Weights;
  float        m_RobustPCAPrecision;
  int Counter;

  vtkLandmarkTransform * m_LandmarkTransform;

private:
  RobustStatisticalShapeModel(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRobustStatisticalShapeModel.txx"
#endif

#endif
