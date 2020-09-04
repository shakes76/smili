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

#ifndef __itkRobustStatisticalShapeModel_txx
#define __itkRobustStatisticalShapeModel_txx

#include "itkRobustStatisticalShapeModel.h"

#include <fstream>
#include <iostream>

#include <vtkFloatArray.h>
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include <vtkMatrix4x4.h>

#include <vtkSmartPointer.h>
#include <vtkLandmarkTransform.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkPolyDataWriter.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkPointData.h>

namespace itk
{

/**
 *
 */
template <class TProfileSamplingPrecisionType>
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::RobustStatisticalShapeModel()
{
  m_PCA = vtkWeightedPCAAnalysisFilter::New();
  m_LandmarkTransform = vtkLandmarkTransform::New();
//  m_PCA->DebugOn();
  m_RobustPCAPrecision = 0.95;
  m_Loaded = false;
  m_PreAligned = false;
  m_ScaleInvariance = true;
  m_MissingMode = true;
  m_RandomInitialisation = false;
  m_Tolerance = 1e-9;
  m_MaxIterations = 3;

  this->m_WeightMode = false;
  this->m_RobustMode = true;
  this->m_StandardMode = false;
}

template <class TProfileSamplingPrecisionType>
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::~RobustStatisticalShapeModel()
{
  if (m_PCA)
    m_PCA->Delete();
  if (m_LandmarkTransform)
    m_LandmarkTransform->Delete();
  ///\todo Check deletion of shapes, original SSM class deletes all shapes itself, this is different for Robust SSM class.
  this->RemoveAllWeights();
  itkDebugMacro(<<"RobustStatisticalShapeModel Destructor");
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GenerateData()
{
  itkDebugMacro("RSSM: GenerateData");
  this->UpdateProcrustesAlignment();
  this->UpdatePCA();

  if (this->m_MeanShape != NULL)
    this->m_MeanShape->Delete();
  this->m_MeanShape = vtkPolyData::New();
  this->m_MeanShape->DeepCopy(this->GetShape(0));
  vtkFloatArray* b = vtkFloatArray::New();
  m_PCA->GetParameterisedShape(b, this->m_MeanShape);
  b->Delete();
  this->m_MeanScale = (double)(this->GetCentroidSize(this->m_MeanShape));

  this->SetValid(true);
}

template <class TProfileSamplingPrecisionType>
bool
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::AddShape(vtkPolyData *shape, vtkFloatArray *weights)
{
  this->SetValid(false);
  // TODO: Add check if shape has same number of points and vertices
  // TODO: Maybe should also check it doesn't already exist !
  //m_Landmarks->AddItem(shape);
  this->m_Landmarks.push_back(shape);
  m_Weights.push_back(weights);
  return true;
}

template <class TProfileSamplingPrecisionType>
bool
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::SetShape(vtkPolyData *shape, vtkFloatArray *weights, unsigned int i)
{
  this->SetValid(false);
  // TODO: Add check if shape has same number of points and vertices
  // TODO: Maybe should also check it doesn't already exist !
  //m_Landmarks->AddItem(shape);
  this->m_Landmarks[i] = shape;
  m_Weights[i] = weights;
  return true;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData*
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetWeightedShape(unsigned int i) throw(itk::ExceptionObject)
{
  if ( i < this->m_Landmarks.size())
  {
    //Get lowest value on mesh
    double range[2];
    m_Weights[i]->GetRange(range);
    double thresholdValue = 1e-3; //clipped all parts above minimum value
    if(range[0] == range[1])
      thresholdValue = range[0]; //no focusing

    vtkPolyData *clippedSurfaceX = vtkPolyData::New();
    this->m_Landmarks[i]->GetPointData()->SetScalars(m_Weights[i]);
    vtkSmartPointer<vtkIdTypeArray> idMap = vtkSmartPointer<vtkIdTypeArray>::New();
    this->ClipSurfaceBasedOnScalars(this->m_Landmarks[i], clippedSurfaceX, idMap, thresholdValue);
    return clippedSurfaceX;
  }
  else
    {
    itkExceptionMacro(<< "Index " << i << " out of range");
    }
}

template <class TProfileSamplingPrecisionType>
vtkFloatArray*
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetWeights(unsigned int i) throw(itk::ExceptionObject)
{
  if ( i < m_Weights.size())
    return m_Weights[i];
  else
  {
    itkExceptionMacro(<< "Index " << i << " out of range");
  }
}

template <class TProfileSamplingPrecisionType>
vtkFloatArray *
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetShapeParameters(vtkPolyData *shape, vtkFloatArray *weights, int bsize)
{
  this->Update();

  vtkPolyData *newShape = vtkPolyData::New(); //overwritten
  newShape->SetPoints(shape->GetPoints());
  vtkFloatArray * b = vtkFloatArray::New();
  if (bsize == -1)
    bsize = this->m_PCA->GetModesRequiredFor(this->GetPrecision());
//  cout << "Using " << bsize << " modes" << endl;
  b->SetNumberOfValues(bsize);
  m_PCA->GetShapeParameters(newShape, weights, b, bsize);

  return b;
}

template <class TProfileSamplingPrecisionType>
vtkFloatArray *
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetShapeParametersInPlace(vtkPolyData *shape, vtkFloatArray *weights, int bsize)
{
  this->Update();

  vtkFloatArray * b = vtkFloatArray::New();
  if (bsize == -1)
    bsize = this->m_PCA->GetModesRequiredFor(this->GetPrecision());
//  cout << "Using " << bsize << " modes" << endl;
  b->SetNumberOfValues(bsize);
  m_PCA->GetShapeParameters(shape, weights, b, bsize);

  return b;
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::RemoveShape(unsigned int id)
{
  itkDebugMacro(<< "Removing Shape " << id);
  unsigned int sz = this->m_Landmarks.size();
  if (id < sz)
  {
    this->m_Landmarks[id]->Delete();
    this->m_Landmarks.erase(this->m_Landmarks.begin() + id);
    this->m_ProcrustesAlignedPoints[id]->Delete();
    this->m_ProcrustesAlignedPoints.erase(this->m_ProcrustesAlignedPoints.begin() + id);
    m_Weights[id]->Delete();
    m_Weights.erase(m_Weights.begin() + id);
  }
  else
    std::cout << "Id out of Range " << id << std::endl;
  if (this->m_Landmarks.size() != (sz-1))
    std::cout << "Failed to delete ?" << std::endl;
  this->SetValid(false);
  this->m_Loaded = false;
  if(this->m_PCA)
    this->m_PCA->Delete();
  this->m_PCA = vtkWeightedPCAAnalysisFilter::New(); //reset
  this->m_PCA->PreLoadedOff();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::RemoveAllWeights()
{
  itkDebugMacro(<< "removing all weights " <<  m_Weights.size());
  int sz = m_Weights.size();
  for (int i = sz-1; i >= 0; i--)
  {
  	//need check if the poly data point to same scalars array and delete only if not
  	//otherwise double deletion of the array when destruction polyDatas in superclass
  	//destructor
    if(m_Weights[i])
      m_Weights[i]->Delete();
    m_Weights.erase(m_Weights.begin() + i);
  }
  if (m_Weights.size() != 0)
    std::cout << "Failed to delete all items ?" << std::endl;
  this->SetValid(false);
}

template <class TProfileSamplingPrecisionType>
float
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetPrecisionRequiredFrom(int modes)
{
  const float resolution = 0.001;

  for(float proportion = resolution; proportion < 1.0; proportion += resolution)
  {
    if(this->GetPCA()->GetModesRequiredFor(proportion) == modes)
      return proportion;
  }

  return 1.0;
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::UpdateProcrustesAlignment() throw(itk::ExceptionObject)
{
  if(m_Loaded)
    return;

  vtkProcrustesAlignmentFilter *procrustesAlign = vtkProcrustesAlignmentFilter::New();

  this->RemoveProcrustesAlignedPoints();
  std::cout << "Computing Weighted Procrustes Alignment" << std::endl;

  // Use rigid body alignment
  if(this->m_PoseType == 1)
    procrustesAlign->GetLandmarkTransform()->SetModeToRigidBody();
  if(this->m_PoseType == 2)
    procrustesAlign->GetLandmarkTransform()->SetModeToSimilarity();
  else if(this->m_PoseType == 3)
    procrustesAlign->GetLandmarkTransform()->SetModeToAffine();
  else
    procrustesAlign->GetLandmarkTransform()->SetModeToRigidBody();

  int sz = this->m_Landmarks.size();
#if VTK_MAJOR_VERSION <= 5
  procrustesAlign->SetNumberOfInputs(sz);
#endif
  procrustesAlign->StartFromCentroidOn();
  for(int i = 0; i < sz; i++)
    {
    vtkPolyData *clipSurface = this->GetWeightedShape(i);
#if VTK_MAJOR_VERSION <= 5
    procrustesAlign->SetInput(i, clipSurface);
#else
    procrustesAlign->AddInputDataObject(i, clipSurface);
#endif
    clipSurface->Delete();
    }
  try
    {
    procrustesAlign->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    itkExceptionMacro(<< "Failed Procrustes Align " << err);
    }

  for(int i = 0; i < sz; i++)
    {
    //Get lowest value on mesh
    double range[2];
    m_Weights[i]->GetRange(range);
    double thresholdValue = 1e-3; //clipped all parts above minimum value
    if(range[0] == range[1])
      thresholdValue = range[0]; //no focusing

    //clip
    this->m_Landmarks[i]->GetPointData()->SetScalars(m_Weights[i]);
    vtkSmartPointer<vtkPolyData> clippedSurfaceX = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkIdTypeArray> idMap = vtkSmartPointer<vtkIdTypeArray>::New();
    this->ClipSurfaceBasedOnScalars(this->m_Landmarks[i], clippedSurfaceX, idMap, thresholdValue);
//    std::cout << "Number of Points: " << clippedSurfaceX->GetNumberOfPoints() << std::endl;

    //compute transform to weight shape from procrustes
    vtkLandmarkTransform * landmarkTransform = vtkLandmarkTransform::New();
    landmarkTransform->SetModeToSimilarity();
    landmarkTransform->SetSourceLandmarks(clippedSurfaceX->GetPoints());
//    landmarkTransform->SetSourceLandmarks(procrustesAlign->GetOutput(i)->GetPoints());
    landmarkTransform->SetTargetLandmarks(vtkPolyData::SafeDownCast(procrustesAlign->GetOutput(i))->GetPoints());
//    landmarkTransform->SetTargetLandmarks(clippedSurfaceX->GetPoints());
    landmarkTransform->Update();
    vtkMatrix4x4 * matrix = landmarkTransform->GetMatrix();

    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->SetMatrix(matrix);

    vtkSmartPointer<vtkTransformPolyDataFilter> transformSurface = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformSurface->SetTransform(transform);
#if VTK_MAJOR_VERSION <= 5
    transformSurface->SetInput(this->m_Landmarks[i]);
#else
    transformSurface->SetInputData(this->m_Landmarks[i]);
#endif
    transformSurface->Update();

    vtkPolyData * aligned = vtkPolyData::New();
    vtkPoints * points = vtkPoints::New();
    points->DeepCopy(transformSurface->GetOutput()->GetPoints());
    aligned->SetPoints(points);
    points->Delete();
    landmarkTransform->Delete();
    this->m_ProcrustesAlignedPoints.push_back(aligned);
    }
  procrustesAlign->Delete();

  itkDebugMacro(<< "Performed Procrustes Alignment with " << sz << " shapes");
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::UpdatePCA() throw(itk::ExceptionObject)
{
  itkDebugMacro(<< "Update PCA");
  int sz = this->m_ProcrustesAlignedPoints.size();

  if (m_Loaded)
    m_PCA->PreLoadedOn();
  else
  {
  #if VTK_MAJOR_VERSION <= 5
    m_PCA->SetNumberOfInputs(sz);
  #endif
    m_PCA->SetConvergenceTolerance(m_Tolerance);
    m_PCA->SetMaximumIterations(m_MaxIterations);
    m_PCA->SetPrecision(m_RobustPCAPrecision);
    m_PCA->SetMissingMode(m_MissingMode);
    if(m_RandomInitialisation)
      m_PCA->RandomGuessOn();
    else
      m_PCA->RandomGuessOff();
  //    m_PCA->DebugOn();
    for (int i = 0; i < sz; i++)
    {
  #if VTK_MAJOR_VERSION <= 5
      m_PCA->SetInput(i, this->m_ProcrustesAlignedPoints[i], m_Weights[i]);
  #else
      m_PCA->SetInput(i, this->m_ProcrustesAlignedPoints[i], m_Weights[i]);
  #endif
    }
  //    cerr << "Weights Set: \n" << m_PCA->GetWeights() << endl;
  }

  try
  {
    m_PCA->Update();
  }
  catch ( itk::ExceptionObject & err )
  {
    itkExceptionMacro(<< "Failed PCA " << err);
  }

//  cerr << "(Weighted-Unweighted) Mean: \n" << m_PCA->GetMeanShape() - m_PCA->GetUnweightedMeanShape() << endl;

  this->m_Modes = m_PCA->GetModesRequiredFor(this->m_Precision);
  itkDebugMacro("Performed PCA Analysis with " << this->m_Modes << " significant modes");
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetParameterisedShape(vtkFloatArray *b)
{
  this->Update();
  //std::cout << "Number of tuples in b is " << b->GetNumberOfTuples() << std::endl;
  // Create a set of pts the same size as the number of landmarks
  vtkPolyData * shape = vtkPolyData::New();
  shape->DeepCopy(this->GetMeanShape());

  m_PCA->GetParameterisedShape(b, shape);

  return shape;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetParameterisedShape(vtkFloatArray *b, vtkFloatArray *weights)
{
  this->Update();
  //std::cout << "Number of tuples in b is " << b->GetNumberOfTuples() << std::endl;
  // Create a set of pts the same size as the number of landmarks
  vtkPolyData * shape = vtkPolyData::New();
  shape->DeepCopy(this->GetMeanShape());

  m_PCA->GetParameterisedShape(b, shape, weights);

  return shape;
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurfaceSimilarityParameters(vtkPolyData * surface, vtkFloatArray *weights, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b)
{
  this->GetSurfaceParametersToMatchShape1(surface, weights, scale, tx, ty, tz, theta, phi, psi, b);
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurfaceSimilarityParameters(vtkPolyData * surface, vtkFloatArray *weights, vtkMatrix4x4 *matrix, vtkFloatArray *b)
{
  itkDebugMacro("GetSurfaceSimilarityParameters with matrix");
  if(m_MissingMode)
    this->GetSurfaceParametersToMatchShape2(surface, weights, matrix, b);
  else
    this->GetSurfaceParametersToMatchShape3(surface, weights, matrix, b);
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurfaceParametersToMatchShape1(vtkPolyData * surface, vtkFloatArray *weights, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b)
{
  // Given a surface find the parameters that minimises the distance between the surface and the generated surfaceArea
  itkDebugMacro("GetSurfaceParametersToMatchShape1");
  bool converged = false;
  int maximumIteration = 50;
  vtkPolyData * generatedSurface = NULL;
  vtkPolyData * oldGeneratedSurface = NULL;
  int i = 0;

  b->FillComponent(0, 0.0);
  theta = 0;
  phi = 0;
  psi = 0;
  Counter = 0;

  while ((converged == false) && (i < maximumIteration))
  {
    //~ std::ostringstream toString;
    //~ toString << Counter;
    //~ std::string iter = "tmp/match_" + toString.str() + "_0.vtk";
    //~ vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
      //~ writer->SetInput(surface);
      //~ writer->SetFileName(iter.c_str());
      //~ writer->Write();

    // Given surface and b -> return pose
    //~ this->GetSurfacePose(surface, scale, tx, ty, tz, theta, phi, psi, b);
    this->GetMissingSurfacePose2(surface, weights, scale, tx, ty, tz, theta, phi, psi, b);

    // Give surface, filter out pose -> return b
    //~ this->GetSurfaceShapeParams(surface, weights, scale, tx, ty, tz, theta, phi, psi, b);
    this->GetMissingSurfaceShapeParams2(surface, weights, scale, tx, ty, tz, theta, phi, psi, b);

    // restrict shape params to +- 3
    double * value2 = new double[1];
    for (int j = 0; j < b->GetNumberOfTuples(); j++)
    {
      b->GetTuple(j, value2);
      if (value2[0] > 3)
        value2[0] = 3;
      else if (value2[0] < -3)
        value2[0] = -3;
      b->SetTuple(j,value2);
    }
    delete [] value2;

    //~ generatedSurface = this->GetSurface(scale, tx, ty, tz, theta, phi, psi, b);
    vtkMatrix4x4 * matrix = m_LandmarkTransform->GetMatrix();
    matrix->Invert(); //From canonical domain to real space
    generatedSurface = this->GetSurface(matrix, b, false);

    //~ iter = "tmp/match_" + toString.str() + "_X.vtk";
    //~ vtkSmartPointer<vtkPolyDataWriter> writer2 = vtkSmartPointer<vtkPolyDataWriter>::New();
      //~ writer2->SetInput(generatedSurface);
      //~ writer2->SetFileName(iter.c_str());
      //~ writer2->Write();

    //Replace missing data with reconstructed ones
    //~ if(hasMissingData)
    //~ {
      //~ for(int j = 0; j < weights->GetNumberOfTuples(); j++)
      //~ {
        //~ if(weights->GetValue(j) == 0)
          //~ surface->GetPoints()->SetPoint(j, generatedSurface->GetPoints()->GetPoint(j));
      //~ }
    //~ }

    if (i > 0)
    {
      converged = this->CheckConvergence(generatedSurface, oldGeneratedSurface, 1e-6);
    }
    if (oldGeneratedSurface != NULL)
      oldGeneratedSurface->Delete();
    oldGeneratedSurface = vtkPolyData::New();
    oldGeneratedSurface->DeepCopy(generatedSurface);
    generatedSurface->Delete();

    i++;
    Counter ++;
  }

  //Get the actual real pose
  //~ this->GetSurfaceShapeParams(surface, weights, scale, tx, ty, tz, theta, phi, psi, b);

  oldGeneratedSurface->Delete();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurfaceParametersToMatchShape2(vtkPolyData * surface, vtkFloatArray *weights, vtkMatrix4x4 *matrix, vtkFloatArray *b)
{
  // Given a surface find the parameters that minimises the distance between the surface and the generated surfaceArea
  itkDebugMacro("GetSurfaceParametersToMatchShape2");
  bool converged = false;
  int maximumIteration = 50;
  vtkPolyData * generatedSurface = NULL;
  vtkPolyData * oldGeneratedSurface = NULL;
  int i = 0;

  b->FillComponent(0, 0.0);
  Counter = 0;

  if(surface->GetNumberOfPoints() == 0)
  {
    itkExceptionMacro(<< "Surface has no points!");
    return;
  }

  while ((converged == false) && (i < maximumIteration))
  {
    surface->Modified();

    // Given surface and b -> return pose
    this->GetMissingSurfacePose(surface, weights, b);

    // Give surface, filter out pose -> return b
    this->GetMissingSurfaceShapeParams(surface, weights, b);

    // restrict shape params to +- 3
    double * value2 = new double[1];
    for (int j = 0; j < b->GetNumberOfTuples(); j++)
    {
      b->GetTuple(j, value2);
      if (value2[0] > 3)
        value2[0] = 3;
      else if (value2[0] < -3)
        value2[0] = -3;
      b->SetTuple(j,value2);
    }
    delete [] value2;

    matrix->DeepCopy(m_LandmarkTransform->GetMatrix());
    matrix->Invert(); //From canonical domain to real space
    generatedSurface = this->GetSurface(matrix, b, false);

    //Replace missing data with reconstructed ones
    for(int j = 0; j < weights->GetNumberOfTuples(); j++)
    {
      if(weights->GetValue(j) == 0)
        surface->GetPoints()->SetPoint(j, generatedSurface->GetPoints()->GetPoint(j));
    }

    if (i > 0)
    {
      converged = this->CheckConvergence(generatedSurface, oldGeneratedSurface, 1e-6);
    }
    if (oldGeneratedSurface != NULL)
      oldGeneratedSurface->Delete();
    oldGeneratedSurface = vtkPolyData::New();
    oldGeneratedSurface->DeepCopy(generatedSurface);
    generatedSurface->Delete();

    i++;
    Counter ++;
  }

  oldGeneratedSurface->Delete();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurfaceParametersToMatchShape3(vtkPolyData * surface, vtkFloatArray *weights, vtkMatrix4x4 *matrix, vtkFloatArray *b)
{
  // Given a surface find the parameters that minimises the distance between the surface and the generated surfaceArea
  itkDebugMacro("GetSurfaceParametersToMatchShape3");
//  std::cerr << "GetSurfaceParametersToMatchShape3" << std::endl;
  bool converged = false;
  int maximumIteration = 50;
  vtkPolyData * generatedSurface = NULL;
  vtkPolyData * oldGeneratedSurface = NULL;
  int i = 0;

  b->FillComponent(0, 0.0);
  Counter = 0;

  if(surface->GetNumberOfPoints() == 0)
  {
    itkExceptionMacro(<< "Surface has no points!");
    return;
  }

  while ((converged == false) && (i < maximumIteration))
  {
    surface->Modified();

//     std::ostringstream toString;
//     toString << Counter;
//     std::string iter = "tmp/match_" + toString.str() + "_0.vtk";
//     vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
//       writer->SetInput(surface);
//       writer->SetFileName(iter.c_str());
//       writer->Write();

    // Given surface and b -> return pose
    this->GetWeightedSurfacePose(surface, weights, b);

    // Give surface, filter out pose -> return b
    this->GetWeightedSurfaceShapeParams(surface, weights, b);

    // restrict shape params to +- 3
    double * value2 = new double[1];
    for (int j = 0; j < b->GetNumberOfTuples(); j++)
    {
      b->GetTuple(j, value2);
      if (value2[0] > 3)
        value2[0] = 3;
      else if (value2[0] < -3)
        value2[0] = -3;
      b->SetTuple(j,value2);
    }
    delete [] value2;

    matrix->DeepCopy(m_LandmarkTransform->GetMatrix());
    matrix->Invert(); //From canonical domain to real space
    generatedSurface = this->GetSurface(matrix, b, false);

//     iter = "tmp/match_" + toString.str() + "_X.vtk";
//     vtkSmartPointer<vtkPolyDataWriter> writer2 = vtkSmartPointer<vtkPolyDataWriter>::New();
//       writer2->SetInput(generatedSurface);
//       writer2->SetFileName(iter.c_str());
//       writer2->Write();

    if (i > 0)
    {
      converged = this->CheckConvergence(generatedSurface, oldGeneratedSurface, 1e-6);
    }
    if (oldGeneratedSurface != NULL)
      oldGeneratedSurface->Delete();
    oldGeneratedSurface = vtkPolyData::New();
    oldGeneratedSurface->DeepCopy(generatedSurface);
    generatedSurface->Delete();

    i++;
    Counter ++;
  }

  oldGeneratedSurface->Delete();
}

template <class TProfileSamplingPrecisionType>
double
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::MahalanobisDistanceBetween(vtkPolyData * poly1, vtkFloatArray *weights1, vtkPolyData * poly2, vtkFloatArray *weights2)
{
  //Setup eigenfeature vectors
  vtkFloatArray *b1 = vtkFloatArray::New();
    b1->SetNumberOfComponents(1);
    b1->SetNumberOfTuples(this->GetNumberOfModes());
    b1->FillComponent(0, 0.0);
  vtkFloatArray *b2 = vtkFloatArray::New();
    b2->SetNumberOfComponents(1);
    b2->SetNumberOfTuples(this->GetNumberOfModes());
    b2->FillComponent(0, 0.0);

  // Now find shape parameters by stripping pose parameters and projecting shapes into model
  vtkSmartPointer<vtkMatrix4x4> matrix1 = vtkSmartPointer<vtkMatrix4x4>::New();
  vtkSmartPointer<vtkMatrix4x4> matrix2 = vtkSmartPointer<vtkMatrix4x4>::New();
  this->GetSurfaceSimilarityParameters(poly1, weights1, matrix1, b1);
  this->GetSurfaceSimilarityParameters(poly2, weights2, matrix2, b2);

//  vtkFloatArray *b1 = this->GetShapeParameters(poly1, this->GetNumberOfModes());
//  vtkFloatArray *b2 = this->GetShapeParameters(poly2, this->GetNumberOfModes());

  //Get eigevalues of PCA
  vtkFloatArray * eigenvalues = m_PCA->GetEvals();

  //Compute Norm
  float norm1 = 0.0, norm2 = 0.0;
  for(int i = 0; i < this->GetNumberOfModes(); i ++)
  {
    b1->SetValue(i, b1->GetValue(i)*sqrt(eigenvalues->GetValue(i)) ); //Undo vtk PCA internal scaling
    b2->SetValue(i, b2->GetValue(i)*sqrt(eigenvalues->GetValue(i)) ); //Undo vtk PCA internal scaling

    norm1 += sqrt(b1->GetValue(i)*b1->GetValue(i));
    norm2 += sqrt(b2->GetValue(i)*b2->GetValue(i));
  }

  double distance = 0.0;
  for(int i = 0; i < this->GetNumberOfModes(); i ++)
    distance += b1->GetValue(i)*b2->GetValue(i)/sqrt(eigenvalues->GetValue(i));
//    distance += b1->GetValue(i)*b2->GetValue(i);
  distance /= norm1*norm2; //Norm

  b1->Delete();
  b2->Delete();

  return -distance;
}

template <class TProfileSamplingPrecisionType>
bool
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::LoadModel(const char * filename) throw(itk::ExceptionObject)
{
  itkDebugMacro(<< "Loading Robust SSM with FileName " << filename);

  typedef int headerType; //long long is 64-bit always
  typedef double writeType;
  typedef int settingsType; //integer settings

  std::ifstream fin(filename, std::ios::binary);
  if (!fin.is_open())
  {
    itkExceptionMacro(<< " Robust SSM3D File " << filename << " does not exist ");
    return false;
  }

  int headerSize = 4;
  headerType *header = new headerType[headerSize];
  fin.read((char*)header, headerSize*sizeof(headerType));

  headerType totalShapes = header[0];
  headerType totalShapesPoints = header[0]*header[1]*header[2];
  headerType totalAlignedPoints = header[0]*header[1]*header[2];
  headerType totalMeanPoints = header[1]*header[2];
  headerType totalCellvalues = header[3];
  headerType totalWeightValues = header[0]*header[2];
  headerType totalEigenvalues = totalShapes;
  headerType totalEigenvectors = totalShapesPoints;
  headerType totalInvEigenvectors = totalShapesPoints;
  headerType totalDataMatrix = totalShapesPoints;

  headerType dataSize = totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues
                 + totalWeightValues + totalEigenvalues + totalEigenvectors + totalInvEigenvectors
                 + totalDataMatrix;
  writeType *data = new writeType[dataSize];

  headerType totalSettings = 5;
  settingsType *settings = new settingsType[totalSettings];

  //std::cout << "Amount of data to read is " << dataSize << std::endl;
  itkDebugMacro(<< "Number of shapes is " << header[0]);
  itkDebugMacro(<< "Number of Dimensions " << header[1]);
  itkDebugMacro(<< "Number of Landmarks " << header[2]);
  itkDebugMacro(<< "Number of Faces " << header[3]/3.0);

  if ((header[0] <=0) || (header[1] <= 0) || (header[2] <= 0))
  {
    delete [] header;
    fin.close();
    itkExceptionMacro(<< "Not a valid ssm " << filename);
    return false;
  }
  this->RemoveAllShapes();

  // Check file size to ensure it's valid using header info
  headerType begin;
  begin = fin.tellg();
  fin.seekg (0, std::ios::end);
  fin.seekg(begin);
  fin.read((char*)data, dataSize*sizeof(writeType));
  fin.read((char*)settings, totalSettings*sizeof(settingsType));
  if (!fin)
  {
    delete [] header;
    delete [] data;
    fin.close();
    itkDebugMacro(<< "Tried reading Full Robust SSM data from file. Trying Compact Robust SSM format instead for " << filename);
    cout << "Tried reading Full Robust SSM data from file. Trying Compact Robust SSM format instead for " << filename << endl;
    return LoadCompactModel(filename);
  }
  m_Loaded = true;
  m_PreAligned = true;

  vtkIdType *faceId = new vtkIdType[3];
  int shapesCount = 0, alignedCount = 0, weightsCount = 0;
  // For each shape get PointSet
  itkDebugMacro(<< "Loading Shapes");
#if VTK_MAJOR_VERSION <= 5
  m_PCA->SetNumberOfInputs(totalShapes);
#endif
  for (int i = 0; i < totalShapes; i++)
  {
    vtkPolyData * shape = vtkPolyData::New();
    // Add points to shape i
    //std::cout << "About to add Points to polyData" << std::endl;
    vtkPoints *points=vtkPoints::New();
    for (int j = 0; j < header[2]; j++)
    {
      // Store each point in array before adding to pointset
      writeType *value = new writeType[header[1]];
      for (int k = 0; k < header[1]; k++)
      {
        value[k] = data[shapesCount];
        shapesCount++;
      }
      points->InsertNextPoint(value);
      delete [] value;
    }
    shape->SetPoints(points);
    points->Delete();

    vtkPolyData * aligned = vtkPolyData::New();
    points=vtkPoints::New();
    for (int j = 0; j < header[2]; j++)
    {
      // Store each point in array before adding to pointset
      writeType *value = new writeType[header[1]];
      for (int k = 0; k < header[1]; k++)
      {
        value[k] = data[totalShapesPoints + alignedCount];
        alignedCount++;
      }
      points->InsertNextPoint(value);
      delete [] value;
    }
    aligned->SetPoints(points);
    points->Delete();

    if(i == 0) //Do only once
    {
      if (m_PCA->GetMeanShape().size() == 0)
      {
        vnl_vector<double> meanVector(totalMeanPoints, 0.0);
        for (size_t j = 0; j < meanVector.size(); j++)
        {
          meanVector[j] = data[totalShapesPoints + totalAlignedPoints + j];
        }

        m_PCA->SetMeanShape(meanVector);
      }
    }

    //std::cout << "About to add Faces to polyData" << std::endl;
    // Add faces for shape i
    vtkCellArray *polys = vtkCellArray::New();
    int cellCount = 0;
    for (int j = 0; j < header[3]/3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        faceId[k] = (vtkIdType) data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + cellCount];
        cellCount ++;
      }
      //std::cout << "Loading (" << faceId[0] << "," << faceId[1] << "," << faceId[2] << ")" << std::endl;
      polys->InsertNextCell(3, faceId);
    }
    shape->SetPolys(polys);
    aligned->SetPolys(polys);

    //Get Weights
    vtkFloatArray *weights = vtkFloatArray::New();
    weights->SetNumberOfTuples(header[2]);
    weights->SetNumberOfComponents(1);
    for (int j = 0; j < weights->GetNumberOfTuples(); j ++)
    {
      weights->SetValue(j, data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + weightsCount]);
      weightsCount ++;
      //~ cout << weights->GetValue(j) << ", ";
    }
    //~ cout << endl;

    polys->Delete();
    //std::cout << "Number of points in shape " << i << " is " << shape->GetNumberOfPoints() << std::endl;
    this->AddShape(shape, weights);
    this->AddAlignedShape(aligned);
  #if VTK_MAJOR_VERSION <= 5
    this->m_PCA->SetInput(i, aligned, weights);
  #else
    this->m_PCA->SetInput(i, aligned, weights);
  #endif
  }

  //Get Eigenvalues and vectors.
  //Get Eigenvalues
  vtkFloatArray *eigenVals = vtkFloatArray::New();
  eigenVals->SetNumberOfTuples(totalEigenvalues);
  eigenVals->SetNumberOfComponents(1);
  for (int j = 0; j < eigenVals->GetNumberOfTuples(); j ++)
  {
    eigenVals->SetValue(j, data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + totalWeightValues + j]);
    //~ cout << eigenVals->GetValue(j) << ", ";
  }
  //~ cout << endl;
  m_PCA->SetEigenValues(eigenVals);
  eigenVals->Delete();
  itkDebugMacro(<< "Loaded Eigenvalues: " << m_PCA->GetEigenValues()->GetNumberOfTuples());

  //Get Eigenvectors
  vnl_matrix<double> vectors(static_cast<double>(totalEigenvectors)/totalShapes, totalShapes, 0.0);
  double *vectorsPtr = vectors.data_block();
//  cerr << "Eigenvectors: " << endl;
  for (int j = 0; j < totalEigenvectors; j ++)
  {
    vectorsPtr[j] = data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + totalWeightValues + totalEigenvalues + j];
//    cerr << vectorsPtr[j] << ", ";
  }
//  cerr << endl;
  m_PCA->SetEigenVectors(vectors);
  itkDebugMacro(<< "Loaded Eigenvectors: " << m_PCA->GetEigenVectors().rows() << "x" << m_PCA->GetEigenVectors().cols());

  //Get Inv Eigenvectors
  vnl_matrix<double> invVectors(totalShapes, static_cast<double>(totalEigenvectors)/totalShapes, 0.0);
  vectorsPtr = invVectors.data_block();
//  cerr << "Eigenvectors: " << endl;
  for (int j = 0; j < totalInvEigenvectors; j ++)
  {
    vectorsPtr[j] = data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + totalWeightValues + totalEigenvalues + totalEigenvectors + j];
//    cerr << vectorsPtr[j] << ", ";
  }
//  cerr << endl;
  m_PCA->SetInverseWeightedEigenVectors(invVectors);
  itkDebugMacro(<< "Loaded Inverse Eigenvectors: " << m_PCA->GetInverseWeightedEigenVectors().rows() << "x" << m_PCA->GetInverseWeightedEigenVectors().cols());

  //Get datamatrix
  vnl_matrix<double> dataMatrix(static_cast<double>(totalDataMatrix)/totalShapes, totalShapes, 0.0);
  vectorsPtr = dataMatrix.data_block();
//  cerr << "Eigenvectors: " << endl;
  for (int j = 0; j < totalDataMatrix; j ++)
  {
    vectorsPtr[j] = data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + totalWeightValues + totalEigenvalues + totalEigenvectors + totalInvEigenvectors + j];
//    cerr << vectorsPtr[j] << ", ";
  }
//  cerr << endl;
  m_PCA->SetDataMatrix(dataMatrix);
  itkDebugMacro(<< "Loaded Data: " << m_PCA->GetDataMatrix().rows() << "x" << m_PCA->GetDataMatrix().cols());

  //Settings
  m_PCA->SetMaximumIterations(settings[0]);
  m_PCA->SetReconstructionIterations(settings[1]);
  m_MaxIterations = m_PCA->GetMaximumIterations();
//  std::cout << "Max Iterations Used for stored model: " << m_PCA->GetMaximumIterations() << std::endl;
//  std::cout << "Iterations Used for stored model: " << m_PCA->GetReconstructionIterations() << std::endl;
  if(settings[2] == 1)
  {
    MissingModeOn();
    m_PCA->MissingModeOn();
    this->m_RobustMode = true;
    this->m_WeightMode = false;
  }
  else
  {
    MissingModeOff();
    m_PCA->MissingModeOff();
    this->m_RobustMode = false;
    this->m_WeightMode = true;
  }
  if(settings[3] == 1)
    PreAlignedOn();
  else
    PreAlignedOff();
  if(settings[4] == 1)
    ScaleInvarianceOn();
  else
    ScaleInvarianceOff();

  //Set the weights in the WPCA object
  UpdatePCA();
  itkDebugMacro(<< "Updated PCA");

  //Update mean polydata
  if (this->m_MeanShape != NULL)
    this->m_MeanShape->Delete();
  this->m_MeanShape = vtkPolyData::New();
  this->m_MeanShape->DeepCopy(this->GetShape(0));
  vtkFloatArray* b = vtkFloatArray::New();
  m_PCA->GetParameterisedShape(b, this->m_MeanShape);
  b->Delete();
  this->m_MeanScale = (double)(this->GetCentroidSize(this->m_MeanShape));
  this->m_StandardMode = false;

  delete [] faceId;
  delete [] data;
  delete [] settings;
  delete [] header;
  itkDebugMacro(<< "Finished Reading SSM Data File");
  return true;
}

template <class TProfileSamplingPrecisionType>
bool
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::LoadCompactModel(const char * filename) throw(itk::ExceptionObject)
{
  itkDebugMacro(<< "Loading Robust SSM with FileName " << filename);

  typedef int headerType; //long long is 64-bit always
  typedef double writeType; 
  typedef int settingsType; //integer settings

  std::ifstream fin(filename, std::ios::binary);
  if (!fin.is_open())
  {
    itkExceptionMacro(<< " Robust SSM3D File " << filename << " does not exist ");
    return false;
  }

  int headerSize = 4;
  headerType *header = new headerType[headerSize];
  fin.read((char*)header, headerSize*sizeof(headerType));

  headerType totalShapes = header[0];
  headerType totalShapesPoints = header[0]*header[1]*header[2];
  headerType totalAlignedPoints = header[0]*header[1]*header[2];
  headerType totalMeanPoints = header[1]*header[2];
  headerType totalCellvalues = header[3];
  headerType totalWeightValues = header[0]*header[2];
  headerType totalEigenvalues = totalShapes;
  headerType totalEigenvectors = totalShapesPoints;

  headerType dataSize = totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues
                 + totalWeightValues + totalEigenvalues + totalEigenvectors;
  writeType *data = new writeType[dataSize];

  headerType totalSettings = 5;
  settingsType *settings = new settingsType[totalSettings];

  //std::cout << "Amount of data to read is " << dataSize << std::endl;
  itkDebugMacro(<< "Number of shapes is " << header[0]);
  itkDebugMacro(<< "Number of Dimensions " << header[1]);
  itkDebugMacro(<< "Number of Landmarks " << header[2]);
  itkDebugMacro(<< "Number of Faces " << header[3]/3.0);

  if ((header[0] <=0) || (header[1] <= 0) || (header[2] <= 0))
  {
    delete [] header;
    fin.close();
    itkExceptionMacro(<< "Not a valid ssm " << filename);
    return false;
  }
  this->RemoveAllShapes();

  // Check file size to ensure it's valid using header info
  headerType begin;
  begin = fin.tellg();
  fin.seekg (0, std::ios::end);
  fin.seekg(begin);
  fin.read((char*)data, dataSize*sizeof(writeType));
  fin.read((char*)settings, totalSettings*sizeof(settingsType));
  if (!fin)
  {
    delete [] header;
    delete [] data;
    fin.close();
    itkExceptionMacro(<< "Failed reading SSM data from file " << filename);
    return false;
  }
  m_Loaded = true;
  m_PreAligned = true;

  vtkIdType *faceId = new vtkIdType[3];
  int shapesCount = 0, alignedCount = 0, weightsCount = 0;
  // For each shape get PointSet
  itkDebugMacro(<< "Loading Shapes");
  m_PCA->SetNumberOfInputs(totalShapes);
  for (int i = 0; i < totalShapes; i++)
  {
    vtkPolyData * shape = vtkPolyData::New();
    // Add points to shape i
    //std::cout << "About to add Points to polyData" << std::endl;
    vtkPoints *points=vtkPoints::New();
    for (int j = 0; j < header[2]; j++)
    {
      // Store each point in array before adding to pointset
      writeType *value = new writeType[header[1]];
      for (int k = 0; k < header[1]; k++)
      {
        value[k] = data[shapesCount];
        shapesCount++;
      }
      points->InsertNextPoint(value);
      delete [] value;
    }
    shape->SetPoints(points);
    points->Delete();

    vtkPolyData * aligned = vtkPolyData::New();
    points=vtkPoints::New();
    for (int j = 0; j < header[2]; j++)
    {
      // Store each point in array before adding to pointset
      writeType *value = new writeType[header[1]];
      for (int k = 0; k < header[1]; k++)
      {
        value[k] = data[totalShapesPoints + alignedCount];
        alignedCount++;
      }
      points->InsertNextPoint(value);
      delete [] value;
    }
    aligned->SetPoints(points);
    points->Delete();

    if(i == 0) //Do only once
    {
      if (m_PCA->GetMeanShape().size() == 0)
      {
        vnl_vector<double> meanVector(totalMeanPoints, 0.0);
        for (size_t j = 0; j < meanVector.size(); j++)
        {
          meanVector[j] = data[totalShapesPoints + totalAlignedPoints + j];
        }

        m_PCA->SetMeanShape(meanVector);
      }
    }

    //std::cout << "About to add Faces to polyData" << std::endl;
    // Add faces for shape i
    vtkCellArray *polys = vtkCellArray::New();
    int cellCount = 0;
    for (int j = 0; j < header[3]/3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        faceId[k] = (vtkIdType) data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + cellCount];
        cellCount ++;
      }
      //std::cout << "Loading (" << faceId[0] << "," << faceId[1] << "," << faceId[2] << ")" << std::endl;
      polys->InsertNextCell(3, faceId);
    }
    shape->SetPolys(polys);
    aligned->SetPolys(polys);

    //Get Weights
    vtkFloatArray *weights = vtkFloatArray::New();
    weights->SetNumberOfTuples(header[2]);
    weights->SetNumberOfComponents(1);
    for (int j = 0; j < weights->GetNumberOfTuples(); j ++)
    {
      weights->SetValue(j, data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + weightsCount]);
      weightsCount ++;
      //~ cout << weights->GetValue(j) << ", ";
    }
    //~ cout << endl;

    polys->Delete();
    //std::cout << "Number of points in shape " << i << " is " << shape->GetNumberOfPoints() << std::endl;
    this->AddShape(shape, weights);
    this->AddAlignedShape(aligned);
  #if VTK_MAJOR_VERSION <= 5
    this->m_PCA->SetInput(i, aligned, weights);
  #else
    this->m_PCA->SetInput(i, aligned, weights);
  #endif
  }

  //Get Eigenvalues and vectors.
  //Get Eigenvalues
  vtkFloatArray *eigenVals = vtkFloatArray::New();
  eigenVals->SetNumberOfTuples(totalEigenvalues);
  eigenVals->SetNumberOfComponents(1);
  for (int j = 0; j < eigenVals->GetNumberOfTuples(); j ++)
  {
    eigenVals->SetValue(j, data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + totalWeightValues + j]);
    //~ cout << eigenVals->GetValue(j) << ", ";
  }
  //~ cout << endl;
  m_PCA->SetEigenValues(eigenVals);
  eigenVals->Delete();
  itkDebugMacro(<< "Loaded Eigenvalues: " << m_PCA->GetEigenValues()->GetNumberOfTuples());

  //Get Eigenvectors
  vnl_matrix<double> vectors(static_cast<double>(totalEigenvectors)/totalShapes, totalShapes, 0.0);
  double *vectorsPtr = vectors.data_block();
//  cerr << "Eigenvectors: " << endl;
  for (int j = 0; j < totalEigenvectors; j ++)
  {
    vectorsPtr[j] = data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + totalWeightValues + totalEigenvalues + j];
//    cerr << vectorsPtr[j] << ", ";
  }
//  cerr << endl;
  m_PCA->SetEigenVectors(vectors);
  itkDebugMacro(<< "Loaded Eigenvectors: " << m_PCA->GetEigenVectors().rows() << "x" << m_PCA->GetEigenVectors().cols());

  //Get Inv Eigenvectors
  m_PCA->SetInverseWeightedEigenVectors(vectors.transpose()); //!< Used for projection, naming convention broken here to avoid branch of code for missing mode
  itkDebugMacro(<< "Loaded Inverse Eigenvectors: " << m_PCA->GetInverseWeightedEigenVectors().rows() << "x" << m_PCA->GetInverseWeightedEigenVectors().cols());

  //Settings
  m_PCA->SetMaximumIterations(settings[0]);
  m_PCA->SetReconstructionIterations(settings[1]);
  m_MaxIterations = m_PCA->GetMaximumIterations();
//  std::cout << "Max Iterations Used for stored model: " << m_PCA->GetMaximumIterations() << std::endl;
//  std::cout << "Iterations Used for stored model: " << m_PCA->GetReconstructionIterations() << std::endl;
  if(settings[2] == 1)
  {
    MissingModeOn();
    m_PCA->MissingModeOn();
    this->m_RobustMode = true;
    this->m_WeightMode = false;
  }
  else
  {
    MissingModeOff();
    m_PCA->MissingModeOff();
    this->m_RobustMode = false;
    this->m_WeightMode = true;
  }
  if(settings[3] == 1)
    PreAlignedOn();
  else
    PreAlignedOff();
  if(settings[4] == 1)
    ScaleInvarianceOn();
  else
    ScaleInvarianceOff();

  //Set the weights in the WPCA object
  UpdatePCA();
  itkDebugMacro(<< "Updated PCA");

  //Update mean polydata
  if (this->m_MeanShape != NULL)
    this->m_MeanShape->Delete();
  this->m_MeanShape = vtkPolyData::New();
  this->m_MeanShape->DeepCopy(this->GetShape(0));
  vtkFloatArray* b = vtkFloatArray::New();
  m_PCA->GetParameterisedShape(b, this->m_MeanShape);
  b->Delete();
  this->m_MeanScale = (double)(this->GetCentroidSize(this->m_MeanShape));
  this->m_StandardMode = false;

  delete [] faceId;
  delete [] data;
  delete [] settings;
  delete [] header;
  itkDebugMacro(<< "Finished Reading SSM Data File");
  return true;
}

template <class TProfileSamplingPrecisionType>
bool
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::SaveModel(const char * filename)
{
  itkDebugMacro(<< "About to Save Robust SSM with filename " << filename);
  /// We want to save the current model to file
  /// File format is header, shape (point) data, aligned (point) data, connectiveness (triangle) data
  typedef int headerType; //vtkIdType is 64-bit always
  typedef double writeType; 
  typedef int settingsType; //integer settings

  headerType sz = this->m_Landmarks.size();
  if (sz > 0)
  {
    std::ofstream fout(filename, std::ios::binary);
    if (!fout.is_open())
    {
      itkExceptionMacro(<< "Cannot create SSM3D File " << filename);
      return false;
    }

    int headerSize = 4;
    headerType * header = new headerType[headerSize];
    header[0] = sz;
    header[1] = 3;
    header[2] = this->GetShape(0)->GetPoints()->GetNumberOfPoints();
    header[3] = 3*this->GetShape(0)->GetPolys()->GetNumberOfCells();

    headerType totalShapes = header[0];
    headerType totalShapesPoints = header[0]*header[1]*header[2];
    headerType totalAlignedPoints = header[0]*header[1]*header[2];
    headerType totalMeanPoints = header[1]*header[2];
    headerType totalCellvalues = header[3];
    headerType totalWeightValues = header[0]*header[2];
    headerType totalEigenvalues = totalShapes;
    headerType totalEigenvectors = totalShapesPoints;
    headerType totalInvEigenvectors = totalShapesPoints;
    headerType totalDataMatrix = totalShapesPoints;

    itkDebugMacro(<< "Header is " << header[0] << " " << header[1] << " " << header[2] << " " << header[3]);
    headerType sizeData = totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues
                          + totalWeightValues + totalEigenvalues + totalEigenvectors + totalInvEigenvectors
                          + totalDataMatrix;
    writeType *writer = new writeType [sizeData];

    headerType totalSettings = 5;
    settingsType *settings = new settingsType [totalSettings];

    ///Extract shape points
    itkDebugMacro(<< "Extracting Shapes.");
    //cerr<< "Extracting Shapes.";
    size_t count = 0;
    writeType value[3];
    for (int i = 0; i < sz; i++)
    {
      vtkPolyData * shape = this->GetShape(i);
      for (int j = 0; j < shape->GetPoints()->GetNumberOfPoints(); j++)
      {
        shape->GetPoints()->GetPoint(j, value);
        writer[count] = value[0];
        count++;
        writer[count] = value[1];
        count++;
        writer[count] = value[2];
        count++;
      }
    }

    ///Extract aligned points
    itkDebugMacro(<< "Extracting Aligned Shapes.");
    for (int i = 0; i < sz; i++)
    {
      vtkPoints * points = this->GetAlignedPoints(i);
      for (int j = 0; j < points->GetNumberOfPoints(); j++)
      {
        points->GetPoint(j, value);
        writer[count] = value[0];
        count++;
        writer[count] = value[1];
        count++;
        writer[count] = value[2];
        count++;
      }
    }

    itkDebugMacro(<< "Extracting Mean Shape.");
    for (size_t j = 0; j < m_PCA->GetMeanShape().size(); j++)
    {
      writer[count] = m_PCA->GetMeanShape()[j];
      count++;
    }

    ///Extract Cells
    itkDebugMacro(<< "Extracting Cells.");
    //cerr<< "Extracting Cells. Coun was " << count;
    this->GetShape(0)->GetPolys()->InitTraversal();
    //std::cout << "Number of faces is " << this->GetShape(0)->GetPolys()->GetNumberOfCells() << std::endl;
    vtkCellArray *cells = this->GetShape(0)->GetPolys();
    vtkTriangle *triangle;
    for (int i=0; i< cells->GetNumberOfCells(); i++)
    {
      vtkIdType fcarray[3];
      triangle=(vtkTriangle *)this->GetShape(0)->GetCell(i);
      fcarray[0] = triangle->GetPointId(0);
      fcarray[1] = triangle->GetPointId(1);
      fcarray[2] = triangle->GetPointId(2);
      //cells->GetNextCell(szv,(vtkIdType *)fcarray);
      //std::cout << i << ": (" << fcarray[0] << "," << fcarray[1] << "," << fcarray[0] << "), ";
      writer[count] = fcarray[0];
      count++;
      writer[count] = fcarray[1];
      count++;
      writer[count] = fcarray[2];
      count++;
    }

    ///Extract Weights
    //cerr<< "Extracting Weights.";
    for (int i = 0; i < sz; i++)
    {
      for (int j = 0; j < m_Weights[i]->GetNumberOfTuples(); j++)
      {
        writer[count] = m_Weights[i]->GetValue(j);
        count++;
      }
    }

    ///Extract Eigenvalues
    //cerr<< "Extracting Eigenvalues.";
    vtkFloatArray *eigenVals = m_PCA->GetEigenValues();
    for (int j = 0; j < eigenVals->GetNumberOfTuples(); j++)
    {
      writer[count] = eigenVals->GetValue(j);
      count++;
    }

    ///Extract Eigenvectors
    //cerr<< "Extracting Eigenvectors.";
    vnl_matrix<double> eigenVecs = m_PCA->GetEigenVectors();
    for (size_t j = 0; j < eigenVecs.rows(); j ++)
    {
      for(size_t k = 0; k < eigenVecs.cols(); k ++)
      {
        writer[count] = eigenVecs(j,k);
        count++;
      }
    }

    ///Extract Inv Eigenvectors
    //cerr<< "Extracting Inv Eigenvectors.";
    vnl_matrix<double> invEigenVecs = m_PCA->GetInverseWeightedEigenVectors();
    for (size_t j = 0; j < invEigenVecs.rows(); j ++)
    {
      for(size_t k = 0; k < invEigenVecs.cols(); k ++)
      {
        writer[count] = invEigenVecs(j,k);
        count++;
      }
    }

    ///Extract Data matrix
    //cerr<< "Extracting Data matrix.";
    vnl_matrix<double> dataMatrix = m_PCA->GetDataMatrix();
    for (size_t j = 0; j < dataMatrix.rows(); j ++)
    {
      for(size_t k = 0; k < dataMatrix.cols(); k ++)
      {
        writer[count] = dataMatrix(j,k);
        count++;
      }
    }

    //Settings
    settings[0] = m_PCA->GetMaximumIterations();
    settings[1] = m_PCA->GetReconstructionIterations();
    if(m_MissingMode)
      settings[2] = 1;
    else
      settings[2] = 0;
    if(m_PreAligned)
      settings[3] = 1;
    else
      settings[3] = 0;
    if(m_ScaleInvariance)
      settings[4] = 1;
    else
      settings[4] = 0;

    // Write data to file
    itkDebugMacro(<< "Writing Data to File.");
    //cerr<< "Writing Data to File.";
    fout.write((char *)(header),headerSize*sizeof(headerType));
    fout.write((char *)(writer),sizeData*sizeof(writeType));
    fout.write((char *)(settings),totalSettings*sizeof(settingsType));
    fout.close();
    delete [] header;
    delete [] writer;
    delete [] settings;
    itkDebugMacro(<< "Finished writing file " << filename);
    cerr<< "Finished writing file " << filename << endl;
  }
  else
  {
    itkExceptionMacro(<< "There is no model to save");
    return false;
  }
  return true;
}

template <class TProfileSamplingPrecisionType>
bool
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::SaveCompactModel(const char * filename)
{
  itkDebugMacro(<< "About to Save Robust SSM (in Compact form) with filename " << filename);
  /// We want to save the current model to file
  /// File format is header, shape (point) data, aligned (point) data, connectiveness (triangle) data
  typedef int headerType; //vtkIdType is 64-bit always
  typedef double writeType; 
  typedef int settingsType; //integer settings

  headerType sz = this->m_Landmarks.size();
  if (sz > 0)
  {
    std::ofstream fout(filename, std::ios::binary);
    if (!fout.is_open())
    {
      itkExceptionMacro(<< "Cannot create SSM3D File " << filename);
      return false;
    }

    int headerSize = 4;
    headerType * header = new headerType[headerSize];
    header[0] = sz;
    header[1] = 3;
    header[2] = this->GetShape(0)->GetPoints()->GetNumberOfPoints();
    header[3] = 3*this->GetShape(0)->GetPolys()->GetNumberOfCells();

    headerType totalShapes = header[0];
    headerType totalShapesPoints = header[0]*header[1]*header[2];
    headerType totalAlignedPoints = header[0]*header[1]*header[2];
    headerType totalMeanPoints = header[1]*header[2];
    headerType totalCellvalues = header[3];
    headerType totalWeightValues = header[0]*header[2];
    headerType totalEigenvalues = totalShapes;
    headerType totalEigenvectors = totalShapesPoints;

    itkDebugMacro(<< "Header is " << header[0] << " " << header[1] << " " << header[2] << " " << header[3]);
    headerType sizeData = totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues
                          + totalWeightValues + totalEigenvalues + totalEigenvectors;
    writeType *writer = new writeType [sizeData];

    headerType totalSettings = 5;
    settingsType *settings = new settingsType [totalSettings];

    ///Extract shape points
    itkDebugMacro(<< "Extracting Shapes.");
    //cerr<< "Extracting Shapes.";
    size_t count = 0;
    writeType value[3];
    for (int i = 0; i < sz; i++)
    {
      vtkPolyData * shape = this->GetShape(i);
      for (int j = 0; j < shape->GetPoints()->GetNumberOfPoints(); j++)
      {
        shape->GetPoints()->GetPoint(j, value);
        writer[count] = value[0];
        count++;
        writer[count] = value[1];
        count++;
        writer[count] = value[2];
        count++;
      }
    }

    ///Extract aligned points
    itkDebugMacro(<< "Extracting Aligned Shapes.");
    for (int i = 0; i < sz; i++)
    {
      vtkPoints * points = this->GetAlignedPoints(i);
      for (int j = 0; j < points->GetNumberOfPoints(); j++)
      {
        points->GetPoint(j, value);
        writer[count] = value[0];
        count++;
        writer[count] = value[1];
        count++;
        writer[count] = value[2];
        count++;
      }
    }

    itkDebugMacro(<< "Extracting Mean Shape.");
    for (size_t j = 0; j < m_PCA->GetMeanShape().size(); j++)
    {
      writer[count] = m_PCA->GetMeanShape()[j];
      count++;
    }

    ///Extract Cells
    itkDebugMacro(<< "Extracting Cells.");
    //cerr<< "Extracting Cells. Coun was " << count;
    this->GetShape(0)->GetPolys()->InitTraversal();
    //std::cout << "Number of faces is " << this->GetShape(0)->GetPolys()->GetNumberOfCells() << std::endl;
    vtkCellArray *cells = this->GetShape(0)->GetPolys();
    vtkTriangle *triangle;
    for (int i=0; i< cells->GetNumberOfCells(); i++)
    {
      vtkIdType fcarray[3];
      triangle=(vtkTriangle *)this->GetShape(0)->GetCell(i);
      fcarray[0] = triangle->GetPointId(0);
      fcarray[1] = triangle->GetPointId(1);
      fcarray[2] = triangle->GetPointId(2);
      //cells->GetNextCell(szv,(vtkIdType *)fcarray);
      //std::cout << i << ": (" << fcarray[0] << "," << fcarray[1] << "," << fcarray[0] << "), ";
      writer[count] = fcarray[0];
      count++;
      writer[count] = fcarray[1];
      count++;
      writer[count] = fcarray[2];
      count++;
    }

    ///Extract Weights
    //cerr<< "Extracting Weights.";
    for (int i = 0; i < sz; i++)
    {
      for (int j = 0; j < m_Weights[i]->GetNumberOfTuples(); j++)
      {
        writer[count] = m_Weights[i]->GetValue(j);
        count++;
      }
    }

    ///Extract Eigenvalues
    //cerr<< "Extracting Eigenvalues.";
    vtkFloatArray *eigenVals = m_PCA->GetEigenValues();
    for (int j = 0; j < eigenVals->GetNumberOfTuples(); j++)
    {
      writer[count] = eigenVals->GetValue(j);
      count++;
    }

    ///Extract Eigenvectors
    //cerr<< "Extracting Eigenvectors.";
    vnl_matrix<double> eigenVecs = m_PCA->GetEigenVectors();
    for (size_t j = 0; j < eigenVecs.rows(); j ++)
    {
      for(size_t k = 0; k < eigenVecs.cols(); k ++)
      {
        writer[count] = eigenVecs(j,k);
        count++;
      }
    }

    //Settings
    settings[0] = m_PCA->GetMaximumIterations();
    settings[1] = m_PCA->GetReconstructionIterations();
    if(m_MissingMode)
      settings[2] = 1;
    else
      settings[2] = 0;
    if(m_PreAligned)
      settings[3] = 1;
    else
      settings[3] = 0;
    if(m_ScaleInvariance)
      settings[4] = 1;
    else
      settings[4] = 0;

    // Write data to file
    itkDebugMacro(<< "Writing Data to File.");
    //cerr<< "Writing Data to File.";
    fout.write((char *)(header),headerSize*sizeof(headerType));
    fout.write((char *)(writer),sizeData*sizeof(writeType));
    fout.write((char *)(settings),totalSettings*sizeof(settingsType));
    fout.close();
    delete [] header;
    delete [] writer;
    delete [] settings;
    itkDebugMacro(<< "Finished writing file " << filename);
    cerr<< "Finished writing file " << filename << endl;
  }
  else
  {
    itkExceptionMacro(<< "There is no model to save");
    return false;
  }
  return true;
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurfacePose2(vtkPolyData * surfaceX, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b)
{
  // Translation is simply the centroid of the surface
  PointType centroidX = this->GetCentroid(surfaceX);
  tx = centroidX[0];
  ty = centroidX[1];
  tz = centroidX[2];
  // The optimal scaling factor is simply derived from the centroid size
  double scaleX = this->GetCentroidSize(surfaceX, centroidX);

  // Generate shape we are aligning with
  vtkPolyData * surfaceY;
  if (b == NULL)
  {
    surfaceY = this->GetMeanShape();
  }
  else
  {
    surfaceY = this->GetParameterisedShape(b);
  }

  PointType centroidY = this->GetCentroid(surfaceY);
  double scaleY = this->GetCentroidSize(surfaceY, centroidY);
  //cout << "GetSurfacePose: ScaleY " << scaleY << endl;
  scale = scaleX/scaleY;
  itkDebugMacro(<< "GetSurfacePose: Scale " << scale);

  vtkLandmarkTransform * landmarkTransform = vtkLandmarkTransform::New();
  landmarkTransform->SetModeToSimilarity();
  landmarkTransform->SetSourceLandmarks(surfaceY->GetPoints());
  landmarkTransform->SetTargetLandmarks(surfaceX->GetPoints());
  landmarkTransform->Update();

  vtkMatrix4x4 * matrix = landmarkTransform->GetMatrix();
  /*
    for(int i = 0; i < 4; i++)
      {
        std::cout << matrix->GetElement(i, 0) << "," << matrix->GetElement(i,1) << "," << matrix->GetElement(i,2) << "," << matrix->GetElement(i,3) << std::endl;
      }*/

  for (int i = 0; i < 3; i++)
  {
    double value = matrix->GetElement(i, 0)/scale;
    matrix->SetElement(i, 0, value);
    value = matrix->GetElement(i, 1)/scale;
    matrix->SetElement(i, 1, value);
    value = matrix->GetElement(i, 2)/scale;
    matrix->SetElement(i, 2, value);
  }
  double value2 = matrix->GetElement(0, 3);
  matrix->SetElement(0, 3, value2 - tx);
  value2 = matrix->GetElement(1, 3);
  matrix->SetElement(1, 3, value2 - ty);
  value2 = matrix->GetElement(2, 3);
  matrix->SetElement(2, 3, value2 - tz);
  for (int i = 0; i < 4; i++)
  {
    std::cout << matrix->GetElement(i, 0) << "," << matrix->GetElement(i,1) << "," << matrix->GetElement(i,2) << "," << matrix->GetElement(i,3) << std::endl;
  }

  PointType orientation2;
  orientation2[0] = 0;
  orientation2[1] = 0;
  orientation2[2] = 0;
  this->GetSurfacePoseFromTransformMatrix(matrix, orientation2);

  theta = orientation2[0];
  phi = orientation2[1];
  psi = orientation2[2];
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurfacePose(vtkPolyData * surfaceX, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b)
{
  // Translation is simply the centroid of the surface
  PointType centroidX = this->GetCentroid(surfaceX);
  tx = centroidX[0];
  ty = centroidX[1];
  tz = centroidX[2];
//  cout << "GetSurfacePose: Translation (" << tx << "," << ty << "," << tz << ")" << endl;

  // The optimal scaling factor is simply derived from the centroid size
  double scaleX = this->GetCentroidSize(surfaceX, centroidX);
//  cout << "GetSurfacePose: ScaleX " << scaleX << endl;

  // Generate shape we are aligning with
  vtkPolyData * surfaceY;
  if (b == NULL)
  {
    surfaceY = this->GetMeanShape();
  }
  else
  {
    surfaceY = this->GetParameterisedShape(b);
  }
  PointType centroidY = this->GetCentroid(surfaceY);
  double scaleY = this->GetCentroidSize(surfaceY, centroidY);
  //cout << "GetSurfacePose: ScaleY " << scaleY << endl;

  scale = scaleX/scaleY;
//  cout << "GetSurfacePose: Scale " << scale << endl;

  // Filter out translation and scale effects
  vtkTransform * transformX = vtkTransform::New();
  transformX->PostMultiply();
  transformX->Identity();
  transformX->Translate(-centroidX[0], -centroidX[1], -centroidX[2]);
  transformX->Scale(1/scale, 1/scale, 1/scale);
  // Centroid Match
  vtkTransformPolyDataFilter *transformSurfaceX = vtkTransformPolyDataFilter::New();
  transformSurfaceX->SetTransform(transformX);
#if VTK_MAJOR_VERSION <= 5
  transformSurfaceX->SetInput(surfaceX);
#else
  transformSurfaceX->SetInputData(surfaceX);
#endif
  transformSurfaceX->Update();

  vtkTransform * transformY = vtkTransform::New();
  transformY->PostMultiply();
  transformY->Identity();
  transformY->Translate(-centroidY[0], -centroidY[1], -centroidY[2]);

  vtkTransformPolyDataFilter *transformSurfaceY = vtkTransformPolyDataFilter::New();
  transformSurfaceY->SetTransform(transformY);
#if VTK_MAJOR_VERSION <= 5
  transformSurfaceY->SetInput(surfaceY);
#else
  transformSurfaceY->SetInputData(surfaceY);
#endif
  transformSurfaceY->Update();

  vtkLandmarkTransform * landmarkTransform = vtkLandmarkTransform::New();
  landmarkTransform->SetModeToRigidBody();
  landmarkTransform->SetSourceLandmarks(transformSurfaceY->GetOutput()->GetPoints());
  landmarkTransform->SetTargetLandmarks(transformSurfaceX->GetOutput()->GetPoints());
  //landmarkTransform->MakeTransform();
  landmarkTransform->Update();
  vtkMatrix4x4 * matrix = landmarkTransform->GetMatrix();

  /*
    for(int i = 0; i < 4; i++)
      {
        std::cout << matrix->GetElement(i, 0) << "," << matrix->GetElement(i,1) << "," << matrix->GetElement(i,2) << "," << matrix->GetElement(i,3) << std::endl;
      }
  */

  PointType orientation2;
  orientation2[0] = 0;
  orientation2[1] = 0;
  orientation2[2] = 0;
  this->GetSurfacePoseFromTransformMatrix(matrix, orientation2);

  theta = orientation2[0];
  phi = orientation2[1];
  psi = orientation2[2];

  surfaceY->Delete();
  transformX->Delete();
  transformSurfaceX->Delete();
  transformY->Delete();
  transformSurfaceY->Delete();
  landmarkTransform->Delete();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetMissingSurfacePose2(vtkPolyData * surfaceX, vtkFloatArray *weights, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b)
{
  //~ std::ostringstream toString;
  //~ toString << Counter;
  //~ std::string iter = "tmp/pose_" + toString.str() + "_O.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer->SetInput(surfaceX);
    //~ writer->SetFileName(iter.c_str());
    //~ writer->Write();

  //Clip surface
  vtkSmartPointer<vtkPolyData> clippedSurfaceX = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkIdTypeArray> idMap = vtkSmartPointer<vtkIdTypeArray>::New();
  this->ClipSurfaceBasedOnScalars(surfaceX, clippedSurfaceX, idMap);

  // Translation is simply the centroid of the surface
  PointType centroidFullX = this->GetCentroid(surfaceX);
  //~ PointType centroidX = this->GetCentroid(clippedSurfaceX);
//  cout << "GetSurfacePose: Translation (" << tx << "," << ty << "," << tz << ")" << endl;

  // The optimal scaling factor is simply derived from the centroid size
  double scaleFullX = this->GetCentroidSize(surfaceX, centroidFullX);
  //~ double scaleX = this->GetCentroidSize(clippedSurfaceX, centroidX);
//  cout << "GetSurfacePose: ScaleX " << scaleX << endl;

  // Generate shape we are aligning with
  vtkPolyData * surfaceY;
  if (b == NULL)
  {
    surfaceY = this->GetMeanShape();
  }
  else
  {
    surfaceY = this->GetParameterisedShape(b);
  }
  surfaceY->GetPointData()->SetScalars(surfaceX->GetPointData()->GetScalars());
  //Clip surface
  vtkSmartPointer<vtkPolyData> clippedSurfaceY = vtkSmartPointer<vtkPolyData>::New();
  this->ClipSurfaceBasedOnScalars(surfaceY, clippedSurfaceY, idMap);

  //Mean shape properties
  PointType centroidMean = this->GetCentroid(this->GetMeanShape());
  double scaleMean = this->GetCentroidSize(this->GetMeanShape(), centroidMean);
  PointType centroidFullY = this->GetCentroid(surfaceY);
  double scaleFullY = this->GetCentroidSize(surfaceY, centroidFullY);
  //~ PointType centroidY = this->GetCentroid(clippedSurfaceY);
  //~ double scaleY = this->GetCentroidSize(clippedSurfaceY, centroidY);
  //cout << "GetSurfacePose: ScaleY " << scaleY << endl;

  tx = centroidFullX[0]-centroidMean[0];
  ty = centroidFullX[1]-centroidMean[1];
  tz = centroidFullX[2]-centroidMean[2];

  scale = scaleFullX/scaleMean;
  double scale2 = scaleFullY/scaleMean;
//  cout << "GetSurfacePose: Scale " << scale << endl;

  // Filter out translation and scale effects
  vtkTransform * transformX = vtkTransform::New();
  transformX->PostMultiply();
  transformX->Identity();
  transformX->Translate(-tx, -ty, -tz);
  transformX->Scale(1/scale, 1/scale, 1/scale);
  // Centroid Match
  vtkTransformPolyDataFilter *transformSurfaceX = vtkTransformPolyDataFilter::New();
  transformSurfaceX->SetTransform(transformX);
#if VTK_MAJOR_VERSION <= 5
  transformSurfaceX->SetInput(clippedSurfaceX);
#else
  transformSurfaceX->SetInputData(clippedSurfaceX);
#endif
  //~ transformSurfaceX->SetInput(surfaceX);
  transformSurfaceX->Update();

  //~ iter = "tmp/pose_" + toString.str() + "_X.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer1 = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer1->SetInput(transformSurfaceX->GetOutput());
    //~ writer1->SetFileName(iter.c_str());
    //~ writer1->Write();

  vtkTransform * transformY = vtkTransform::New();
  transformY->PostMultiply();
  transformY->Identity();
  transformY->Translate(-centroidFullY[0]+centroidMean[0], -centroidFullY[1]+centroidMean[1], -centroidFullY[2]+centroidMean[2]);
  transformX->Scale(1/scale2, 1/scale2, 1/scale2);

  vtkTransformPolyDataFilter *transformSurfaceY = vtkTransformPolyDataFilter::New();
  transformSurfaceY->SetTransform(transformY);
#if VTK_MAJOR_VERSION <= 5
  transformSurfaceY->SetInput(clippedSurfaceY);
#else
  transformSurfaceY->SetInputData(clippedSurfaceY);
#endif
  //~ transformSurfaceY->SetInput(surfaceY);
  transformSurfaceY->Update();

  //~ iter = "tmp/pose_" + toString.str() + "_Y.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer2 = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer2->SetInput(transformSurfaceY->GetOutput());
    //~ writer2->SetFileName(iter.c_str());
    //~ writer2->Write();

  vtkLandmarkTransform * landmarkTransform = vtkLandmarkTransform::New();
  landmarkTransform->SetModeToRigidBody();
  landmarkTransform->SetSourceLandmarks(transformSurfaceY->GetOutput()->GetPoints());
  landmarkTransform->SetTargetLandmarks(transformSurfaceX->GetOutput()->GetPoints());
  //landmarkTransform->MakeTransform();
  landmarkTransform->Update();
  vtkMatrix4x4 * matrix = landmarkTransform->GetMatrix();

  PointType orientation2;
  orientation2[0] = 0;
  orientation2[1] = 0;
  orientation2[2] = 0;
  this->GetSurfacePoseFromTransformMatrix(matrix, orientation2);

  theta = orientation2[0];
  phi = orientation2[1];
  psi = orientation2[2];

  surfaceY->Delete();
  transformX->Delete();
  transformSurfaceX->Delete();
  transformY->Delete();
  transformSurfaceY->Delete();
  landmarkTransform->Delete();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetWeightedSurfacePose(vtkPolyData * surfaceX, vtkFloatArray *weights, vtkFloatArray *b)
{
  itkDebugMacro("GetWeightedSurfacePose");
  // Generate shape we are aligning with
  vtkPolyData * surfaceY;
  if (b == NULL)
  {
    surfaceY = this->GetMeanShape();
  }
  else
  {
    surfaceY = this->GetParameterisedShape(b, weights);
  }
  surfaceY->GetPointData()->SetScalars(weights);

  double range[2];
  weights->GetRange(range);
  double thresholdValue;
  if(vnl_math_abs(range[1]-range[0])<1e-6) // uniform weights
  {
    thresholdValue = range[0] - 1e-6; // include the whole shape
  }
  else
  {
    thresholdValue = range[0] + 1e-6; // clipp all parts above minimum value
  }

  ///\todo Intelligent clipping hereb ased on the weights, currently xlipps the minimum value only

  ///Clip to the main (highly) weighted area(s)
  vtkSmartPointer<vtkPolyData> clippedSurfaceX = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkIdTypeArray> idMapX = vtkSmartPointer<vtkIdTypeArray>::New();
  this->ClipSurfaceBasedOnScalars(surfaceX, clippedSurfaceX, idMapX, thresholdValue);
  vtkSmartPointer<vtkPolyData> clippedSurfaceY = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkIdTypeArray> idMapY = vtkSmartPointer<vtkIdTypeArray>::New();
  this->ClipSurfaceBasedOnScalars(surfaceY, clippedSurfaceY, idMapY, thresholdValue);

//  std::ostringstream toString;
//  toString << Counter;
//  std::string iter = "tmp/pose_" + toString.str() + "_clipX.vtk";
//  vtkSmartPointer<vtkPolyDataWriter> writer1 = vtkSmartPointer<vtkPolyDataWriter>::New();
//   writer1->SetInput(clippedSurfaceX);
//   writer1->SetFileName(iter.c_str());
//   writer1->Write();
//
//  iter = "tmp/pose_" + toString.str() + "_clipY.vtk";
//  vtkSmartPointer<vtkPolyDataWriter> writer2 = vtkSmartPointer<vtkPolyDataWriter>::New();
//   writer2->SetInput(clippedSurfaceY);
//   writer2->SetFileName(iter.c_str());
//   writer2->Write();

  m_LandmarkTransform->SetSourceLandmarks(clippedSurfaceX->GetPoints());
  m_LandmarkTransform->SetTargetLandmarks(clippedSurfaceY->GetPoints());
  if(m_ScaleInvariance)
    m_LandmarkTransform->SetModeToSimilarity();
  else
    m_LandmarkTransform->SetModeToRigidBody();
  m_LandmarkTransform->Update();

  surfaceY->Delete();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetMissingSurfacePose(vtkPolyData * surfaceX, vtkFloatArray *weights, vtkFloatArray *b)
{
  itkDebugMacro("GetMissingSurfacePose");
  //Clip original surface
  vtkSmartPointer<vtkPolyData> clippedSurfaceX = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkIdTypeArray> idMap = vtkSmartPointer<vtkIdTypeArray>::New();
  this->ClipSurfaceBasedOnScalars(surfaceX, clippedSurfaceX, idMap);

  // Generate shape we are aligning with
  vtkPolyData * surfaceY;
  if (b == NULL)
  {
    surfaceY = this->GetMeanShape();
  }
  else
  {
    surfaceY = this->GetParameterisedShape(b);
  }
  surfaceY->GetPointData()->SetScalars(surfaceX->GetPointData()->GetScalars());
  //Clip parameterised surface
  vtkSmartPointer<vtkPolyData> clippedSurfaceY = vtkSmartPointer<vtkPolyData>::New();
  this->ClipSurfaceBasedOnScalars(surfaceY, clippedSurfaceY, idMap);

  //~ std::ostringstream toString;
  //~ toString << Counter;
  //~ std::string iter = "tmp/pose_" + toString.str() + "_clipX.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer1 = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer1->SetInput(clippedSurfaceX);
    //~ writer1->SetFileName(iter.c_str());
    //~ writer1->Write();
  //~
  //~ iter = "tmp/pose_" + toString.str() + "_clipY.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer2 = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer2->SetInput(clippedSurfaceY);
    //~ writer2->SetFileName(iter.c_str());
    //~ writer2->Write();

  m_LandmarkTransform->SetSourceLandmarks(clippedSurfaceX->GetPoints());
  m_LandmarkTransform->SetTargetLandmarks(clippedSurfaceY->GetPoints());
  if(m_ScaleInvariance)
    m_LandmarkTransform->SetModeToSimilarity();
  else
    m_LandmarkTransform->SetModeToRigidBody();
  m_LandmarkTransform->Update();

  //Transform
//  for (int k = 0; k < clippedSurfaceX->GetNumberOfPoints(); k ++)
//  {
//    double outPoint[3];
//
//    m_LandmarkTransform->InternalTransformPoint( clippedSurfaceX->GetPoints()->GetPoint(k), outPoint );
//
//    clippedSurfaceX->GetPoints()->SetPoint(k, outPoint);
//  }

  //~ iter = "tmp/pose_" + toString.str() + "_clipX_transformed.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer3 = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer3->SetInput(clippedSurfaceX);
    //~ writer3->SetFileName(iter.c_str());
    //~ writer3->Write();

  surfaceY->Delete();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurfaceShapeParams(vtkPolyData * surface, vtkFloatArray *weights, double scale, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b)
{
  // Remove pose effects -
  vtkTransform *transform = vtkTransform::New();
  /*transform->PostMultiply();
  transform->Identity();
  transform->Translate(tx, ty, tz);
  transform->Scale(scale, scale, scale);
  transform->RotateZ(theta);
  transform->RotateX(phi);
  transform->RotateZ(psi);
  transform->Inverse();*/
  // Reverse transforms

  // Reverse transforms
  // Place
  vtkMatrix4x4 * matrix = this->GetSimilarityMatrix(scale, tx, ty, tz, theta, phi, psi, 0, 0, 0);
  matrix->Invert();

  transform->SetMatrix(matrix);

  //vtkMatrix4x4 * matrix = transform->GetMatrix();
  //cout << "Use Matrix " << endl;
  //cout << matrix->GetElement(0,0) << "," << matrix->GetElement(0,1) << "," << matrix->GetElement(0,2) << endl;
  //cout << matrix->GetElement(1,0) << "," << matrix->GetElement(1,1) << "," << matrix->GetElement(1,2) << endl;
  //cout << matrix->GetElement(2,0) << "," << matrix->GetElement(2,1) << "," << matrix->GetElement(2,2) << endl;

  vtkTransformPolyDataFilter *transformSurface = vtkTransformPolyDataFilter::New();
  transformSurface->SetTransform(transform);
#if VTK_MAJOR_VERSION <= 5
  transformSurface->SetInput(surface);
#else
  transformSurface->SetInputData(surface);
#endif
  transformSurface->Update();

  //~ std::string iter = "GetShapeParams.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer->SetInput(transformSurface->GetOutput());
    //~ writer->SetFileName(iter.c_str());
    //~ writer->Write();

//  PointType centroid = this->GetCentroid(transformSurface->GetOutput());
//  cout << "GetSurfaceShapeParams: Translation (" << centroid[0] << "," << centroid[1] << "," << centroid[2] << ")" << endl;
//  cout << "GetSurfaceShapeParams: Scale " << this->GetCentroidSize(transformSurface->GetOutput()) << endl;

  // Now find shape parameters
  vtkFloatArray * bNew = this->GetShapeParametersInPlace(transformSurface->GetOutput(), weights, b->GetNumberOfTuples());
  b->DeepCopy(bNew);

  bNew->Delete();
  transform->Delete();
  transformSurface->Delete();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetMissingSurfaceShapeParams2(vtkPolyData * surface, vtkFloatArray *weights, double scale, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b)
{
  // Remove pose effects -
  vtkTransform *transform = vtkTransform::New();
  // Reverse transforms
  // Place
  vtkMatrix4x4 * matrix = this->GetSimilarityMatrix(scale, tx, ty, tz, theta, phi, psi, 0, 0, 0);
  matrix->Invert();

  transform->SetMatrix(matrix);

  //~ vtkMatrix4x4 * matrix = transform->GetMatrix();
  //~ cout << "Use Matrix " << endl;
  //~ cout << matrix->GetElement(0,0) << "," << matrix->GetElement(0,1) << "," << matrix->GetElement(0,2) << endl;
  //~ cout << matrix->GetElement(1,0) << "," << matrix->GetElement(1,1) << "," << matrix->GetElement(1,2) << endl;
  //~ cout << matrix->GetElement(2,0) << "," << matrix->GetElement(2,1) << "," << matrix->GetElement(2,2) << endl;

  //Clip surface
  vtkSmartPointer<vtkPolyData> clippedSurface = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkIdTypeArray> idMap = vtkSmartPointer<vtkIdTypeArray>::New();
  this->ClipSurfaceBasedOnScalars(surface, clippedSurface, idMap);

  vtkTransformPolyDataFilter *transformSurface = vtkTransformPolyDataFilter::New();
  transformSurface->SetTransform(transform);
#if VTK_MAJOR_VERSION <= 5
  transformSurface->SetInput(clippedSurface);
#else
  transformSurface->SetInputData(clippedSurface);
#endif
  transformSurface->Update();

  //Restore correspondence
  vtkSmartPointer<vtkPolyData> clippedRemappedSurface = vtkSmartPointer<vtkPolyData>::New();
  this->RestoreClippedSurfaceCorrespondence(transformSurface->GetOutput(), idMap, clippedRemappedSurface);
  clippedRemappedSurface->SetPolys(surface->GetPolys());

  //~ std::ostringstream toString;
  //~ toString << Counter;
  //~ std::string iter = "tmp/param_" + toString.str() + "_X.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer->SetInput(clippedRemappedSurface);
    //~ writer->SetFileName(iter.c_str());
    //~ writer->Write();

  //~ PointType centroid = this->GetCentroid(transformSurface->GetOutput());
  //~ cout << "GetSurfaceShapeParams: Translation (" << centroid[0] << "," << centroid[1] << "," << centroid[2] << ")" << endl;
  //~ cout << "GetSurfaceShapeParams: Scale " << this->GetCentroidSize(transformSurface->GetOutput()) << endl;

  // Now find shape parameters
  vtkFloatArray * bNew = this->GetShapeParametersInPlace(clippedRemappedSurface, weights, b->GetNumberOfTuples());
  b->DeepCopy(bNew);

  //~ iter = "tmp/param_" + toString.str() + "_Y.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer2 = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer2->SetInput(clippedRemappedSurface);
    //~ writer2->SetFileName(iter.c_str());
    //~ writer2->Write();

  bNew->Delete();
  transform->Delete();
  transformSurface->Delete();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetWeightedSurfaceShapeParams(vtkPolyData * surface, vtkFloatArray *weights, vtkFloatArray *b)
{
  itkDebugMacro("GetWeightedSurfaceShapeParams");
  // Reverse transforms
  vtkMatrix4x4 * matrix = m_LandmarkTransform->GetMatrix();
  //~ matrix->Invert();

  // Remove pose effects -
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->SetMatrix(matrix);

  //~ cout << "Use Matrix " << endl;
  //~ cout << matrix->GetElement(0,0) << "," << matrix->GetElement(0,1) << "," << matrix->GetElement(0,2) << endl;
  //~ cout << matrix->GetElement(1,0) << "," << matrix->GetElement(1,1) << "," << matrix->GetElement(1,2) << endl;
  //~ cout << matrix->GetElement(2,0) << "," << matrix->GetElement(2,1) << "," << matrix->GetElement(2,2) << endl;

  vtkSmartPointer<vtkTransformPolyDataFilter> transformSurface = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformSurface->SetTransform(transform);
#if VTK_MAJOR_VERSION <= 5
  transformSurface->SetInput(surface);
#else
  transformSurface->SetInputData(surface);
#endif
  transformSurface->Update();

  // Now find shape parameters
  vtkFloatArray * bNew = this->GetShapeParametersInPlace(transformSurface->GetOutput(), weights, b->GetNumberOfTuples());
  b->DeepCopy(bNew);

  bNew->Delete();
}

template <class TProfileSamplingPrecisionType>
void
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetMissingSurfaceShapeParams(vtkPolyData * surface, vtkFloatArray *weights, vtkFloatArray *b)
{
  itkDebugMacro("GetMissingSurfaceShapeParams");
  // Reverse transforms
  vtkMatrix4x4 * matrix = m_LandmarkTransform->GetMatrix();
  //~ matrix->Invert();

  // Remove pose effects -
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->SetMatrix(matrix);

  //~ vtkMatrix4x4 * matrix = transform->GetMatrix();
  //~ cout << "Use Matrix " << endl;
  //~ cout << matrix->GetElement(0,0) << "," << matrix->GetElement(0,1) << "," << matrix->GetElement(0,2) << endl;
  //~ cout << matrix->GetElement(1,0) << "," << matrix->GetElement(1,1) << "," << matrix->GetElement(1,2) << endl;
  //~ cout << matrix->GetElement(2,0) << "," << matrix->GetElement(2,1) << "," << matrix->GetElement(2,2) << endl;

  //~ std::ostringstream toString;
  //~ toString << Counter;
  //~ std::string iter = "tmp/param_" + toString.str() + "_0.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer->SetInput(surface);
    //~ writer->SetFileName(iter.c_str());
    //~ writer->Write();

  //Clip surface
  vtkSmartPointer<vtkPolyData> clippedSurface = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkIdTypeArray> idMap = vtkSmartPointer<vtkIdTypeArray>::New();
  this->ClipSurfaceBasedOnScalars(surface, clippedSurface, idMap);

  //~ iter = "tmp/param_" + toString.str() + "_X.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer1 = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer1->SetInput(clippedSurface);
    //~ writer1->SetFileName(iter.c_str());
    //~ writer1->Write();

  vtkSmartPointer<vtkTransformPolyDataFilter> transformSurface = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformSurface->SetTransform(transform);
#if VTK_MAJOR_VERSION <= 5
  transformSurface->SetInput(clippedSurface);
#else
  transformSurface->SetInputData(clippedSurface);
#endif
  transformSurface->Update();

  //Restore correspondence
  vtkSmartPointer<vtkPolyData> clippedRemappedSurface = vtkSmartPointer<vtkPolyData>::New();
  this->RestoreClippedSurfaceCorrespondence(transformSurface->GetOutput(), idMap, clippedRemappedSurface);
  clippedRemappedSurface->SetPolys(surface->GetPolys());

  // Now find shape parameters
  vtkFloatArray * bNew = this->GetShapeParametersInPlace(clippedRemappedSurface, weights, b->GetNumberOfTuples());
  b->DeepCopy(bNew);

  //~ iter = "tmp/param_" + toString.str() + "_Y.vtk";
  //~ vtkSmartPointer<vtkPolyDataWriter> writer2 = vtkSmartPointer<vtkPolyDataWriter>::New();
    //~ writer2->SetInput(clippedRemappedSurface);
    //~ writer2->SetFileName(iter.c_str());
    //~ writer2->Write();

  bNew->Delete();
}

template <class TProfileSamplingPrecisionType>
double
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetWeightedCentroidSize(vtkPolyData * surface, vtkFloatArray *weights, PointType centroid)
{
  // TODO: Using true centroid size... but doesn't match definition of scale from Horn so using here.
  PointType pt;
  pt[0] = 0;
  pt[1] = 0;
  pt[2] = 0;
  double centroidSize = 0;
  for (int i = 0; i < surface->GetPoints()->GetNumberOfPoints(); i++)
  {
    if(weights->GetValue(i) != 0.0)
    {
      double *p = surface->GetPoints()->GetPoint(i);
      pt[0] = p[0] - centroid[0];
      pt[1] = p[1] - centroid[1];
      pt[2] = p[2] - centroid[2];
      centroidSize += pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2];
    }
  }
  return sqrt(centroidSize);
}

template <class TProfileSamplingPrecisionType>
typename RobustStatisticalShapeModel<TProfileSamplingPrecisionType>::PointType
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetWeightedCentroid(vtkPolyData * surface, vtkFloatArray *weights)
{
  // Centroid
  PointType centroid;
  centroid[0] = 0;
  centroid[1] = 0;
  centroid[2] = 0;
  int numberPoints = surface->GetPoints()->GetNumberOfPoints();

  // Calculate center of shape
  int counter = 0;
  for (int i = 0; i < numberPoints; i++)
  {
    double *point = surface->GetPoints()->GetPoint(i);
    if(weights->GetValue(i) != 0.0)
    {
      centroid[0] += point[0];
      centroid[1] += point[1];
      centroid[2] += point[2];
      counter ++;
    }
  }
  centroid[0] /= (float)counter;
  centroid[1] /= (float)counter;
  centroid[2] /= (float)counter;

  return centroid;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurface(vtkMatrix4x4 * matrixTransform, vtkFloatArray *b, bool centroid)
{
  vtkPolyData * poly = vtkPolyData::New();

  vtkPolyData * surface;
  if (b == NULL)
    surface = this->GetMeanShape();
  else
    surface = this->GetParameterisedShape(b);

  vtkTransform * transform = vtkTransform::New();
  if (centroid == true)
  {
    PointType centroidPoint = this->GetCentroid(surface);
    transform->PostMultiply();
    transform->Translate(-centroidPoint[0], -centroidPoint[1], -centroidPoint[2]);
    transform->Concatenate(matrixTransform);
  }
  else
  {
    transform->SetMatrix(matrixTransform);
  }
  vtkTransformPolyDataFilter *transformSurface = vtkTransformPolyDataFilter::New();
  transformSurface->SetTransform(transform);
#if VTK_MAJOR_VERSION <= 5
  transformSurface->SetInput(surface);
#else
  transformSurface->SetInputData(surface);
#endif
  transformSurface->Update();

  poly->DeepCopy(transformSurface->GetOutput());
  transformSurface->Delete();
  transform->Delete();
  surface->Delete();

  return poly;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurface(vnl_matrix<double> pose, vtkFloatArray *b, bool centroid)
{
  vtkMatrix4x4 * matrixTransform = vtkMatrix4x4::New();
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      matrixTransform->SetElement(i,j, pose(i,j));
    }
  }

  vtkPolyData * surface = this->GetSurface(matrixTransform, b, centroid);
  matrixTransform->Delete();
  return surface;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurface(double scale, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b)
{
  PointType t;
  t[0] = tx;
  t[1] = ty;
  t[2] = tz;
  PointType orientation;
  orientation[0] = theta;
  orientation[1] = phi;
  orientation[2] = psi;
  vtkPolyData * surface = this->GetSurface(scale, t, orientation, b);
  return surface;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
RobustStatisticalShapeModel<TProfileSamplingPrecisionType>
::GetSurface(double s, PointType t, PointType orientation, vtkFloatArray *b)
{
  vtkMatrix4x4 * matrix = this->GetSimilarityMatrix(s, t[0], t[1], t[2], orientation[0], orientation[1], orientation[2]);
  vtkPolyData * surface =  this->GetSurface(matrix, b, true);
  matrix->Delete();
  return surface;
}

}

#endif
