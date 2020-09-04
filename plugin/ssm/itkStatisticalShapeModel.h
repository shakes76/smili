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

#ifndef __itkStatisticalShapeModel_h
#define __itkStatisticalShapeModel_h

// vtk headers
#include <vtkPCAAnalysisFilter.h>

#include "itkStatisticalShapeModelBase.h"

namespace itk
{

/** 	\class StatisticalShapeModel
	\brief Construct, save and open Statistical Shape Models (SSMs) of Cootes et al. (1995).

	Cootes, T. F.; Taylor, C. J.; Cooper, D. H. & Graham, J.,
	Active Shape Models-Their Training and Application,
	Computer Vision and Image Understanding, 1995, 61, 38 - 59.

	The shapes passed (via AddShape()) are aligned using vtkProcrustesAlignmentFilter and then the PCA
	(via vtkPCAAnalysisFilter) is computed. Shapes can then be parameterised using the principal components of the SSM.

	An example application includes computing the model thickness of cartilage tissue above the bone.

	Declaration Example:
	\code
	typedef itk::StatisticalShapeModel<double> ShapeModelType;

	ShapeModelType::Pointer m_SSM; //!< The Statistical Shape Model
	m_SSM = ShapeModelType::New();
	\endcode

	Usage Example 1:
	\code
	vtkPolyDataCollection* meshes = loadMeshes(); //Load some of your meshes here, implement yourself

	size_t n = meshes->GetNumberOfItems();
	meshes->InitTraversal();
	for(size_t j = 0; j < n; j ++)
		m_SSM->AddShape(meshes->GetNextItem()); //!< Set model to shape filter

	m_SSM->Update();
	//Dont destroy meshes collection, SSM class handles that
	\endcode

	Usage Example 2:
	\code
	m_SSM->LoadModel(filename.toStdString().c_str());
	m_SSM->Update();
	\endcode
 */

template <class TProfileSamplingPrecisionType=double>
class MILX_EXPORT StatisticalShapeModel : public StatisticalShapeModelBase<TProfileSamplingPrecisionType>
{
public:
  typedef StatisticalShapeModel                                         Self;
  typedef StatisticalShapeModelBase<TProfileSamplingPrecisionType>      Superclass;
  typedef SmartPointer<Self>                                            Pointer;
  typedef SmartPointer<const Self>                                      ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(StatisticalShapeModel, StatisticalShapeModelBase);

  /**
   * Load Model from file
   *  - Destroys current model
   */
  virtual bool LoadModel(const char *) throw(itk::ExceptionObject);

  /**
   * Save Model as file
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

  virtual int GetNumberOfModes()
    {
    this->Update();
    return this->GetPCA()->GetModesRequiredFor(this->m_Precision);
    }
  virtual float GetPrecisionRequiredFrom(int modes);

  /**
   * Calculates the Mahalanobis distance between two surfaces within the shape model
   *  - Assumes polyData has same number of points
   *  - Assumes the surfaces have been pre co-registered.
   */
  double MahalanobisDistanceBetween(vtkPolyData * poly1, vtkPolyData * poly2, const bool includeScaling = false);

  /**
   * Generate Shape using b
   */
  virtual vtkPolyData * GetParameterisedShape(vtkFloatArray *b);
  /**
   * Generate Parameters that best describe shape
   */
  vtkFloatArray * GetShapeParameters(vtkPolyData *shape, int bsize = -1);

  /**
    * Allow access to internal structures
    *  - Basically only needed for debugging and visualization (ie should be hidden or removed later)
    */
  inline vtkPCAAnalysisFilter* GetPCA(){ return m_PCA; }

protected:
  StatisticalShapeModel();
  virtual ~StatisticalShapeModel();

  virtual void GenerateData();
  virtual void UpdatePCA() throw(itk::ExceptionObject);

  vtkPCAAnalysisFilter* m_PCA; //!< Computes the actual shape model using PCA

private:
  StatisticalShapeModel(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkStatisticalShapeModel.txx"
#endif

#endif
