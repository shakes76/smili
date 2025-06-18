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

#ifndef __itkStatisticalShapeModelBase_txx
#define __itkStatisticalShapeModelBase_txx

#include "itkStatisticalShapeModelBase.h"

#include <fstream>
#include <iostream>

#include <vtkFloatArray.h>
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include <vtkMatrix4x4.h>

#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkProcrustesAlignmentFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkPolyDataWriter.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkMultiBlockDataSet.h>

namespace itk
{

/**
 *
 */
template <class TProfileSamplingPrecisionType>
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::StatisticalShapeModelBase()
{
  this->m_PoseType = 2;
  m_Precision = 0.9;
  m_MeanScale = 0;
  m_MeanShape = NULL;
  m_Valid = false;

  m_WeightMode = false;
  m_RobustMode = false;
  m_StandardMode = true;
}

template <class TProfileSamplingPrecisionType>
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::~StatisticalShapeModelBase()
{
  if(m_MeanShape != NULL)
    m_MeanShape->Delete();
  this->RemoveAllShapes();
  this->RemoveProcrustesAlignedPoints();
  //m_Landmarks->Delete();
  itkDebugMacro(<<"StatisticalShapeModelBaseBase Destructor");
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::Update()
{
  if(this->GetValid() == false)
    this->GenerateData();
}

template <class TProfileSamplingPrecisionType>
vtkPolyData*
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetMeanShape()
{
  this->Update();
  return m_MeanShape;
}

template <class TProfileSamplingPrecisionType>
bool
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::AddShape(vtkPolyData *shape)
{
  this->SetValid(false);
  // TODO: Add check if shape has same number of points and vertices
  // TODO: Maybe should also check it doesn't already exist !
  //m_Landmarks->AddItem(shape);
  m_Landmarks.push_back(shape);
  return true;
}

template <class TProfileSamplingPrecisionType>
bool
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::SetShape(vtkPolyData *shape, unsigned int i)
{
  this->SetValid(false);
  // TODO: Add check if shape has same number of points and vertices
  // TODO: Maybe should also check it doesn't already exist !
  //m_Landmarks->AddItem(shape);
  m_Landmarks[i] = shape;
  return true;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData*
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetShape(unsigned int i) throw(itk::ExceptionObject)
{
  if ( i < m_Landmarks.size())
    return m_Landmarks[i];
  else
    {
    itkExceptionMacro(<< "Index " << i << " out of range");
    }
}

template <class TProfileSamplingPrecisionType>
bool
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::AddAlignedShape(vtkPolyData *aligned)
{
  this->SetValid(false);
  // TODO: Add check if shape has same number of points and vertices
  // TODO: Maybe should also check it doesn't already exist !
  m_ProcrustesAlignedPoints.push_back(aligned);
  return true;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData*
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetProcrustesAlignedSurface(int index)
{
  this->Update();
  vtkPolyData * shape = this->GetMeanShape();

  vtkPolyData * polyData = vtkPolyData::New();
  vtkCellArray *cells = vtkCellArray::New();

  vtkIdType *faceId = new vtkIdType[3];
  vtkTriangle *triangle;
  shape->GetPolys()->InitTraversal();
  int numberOfCells = shape->GetPolys()->GetNumberOfCells();
  for (int j = 0; j < numberOfCells; j++)
    {
    triangle=(vtkTriangle *)shape->GetCell(j);
    faceId[0] = triangle->GetPointId(0);
    faceId[1] = triangle->GetPointId(1);
    faceId[2] = triangle->GetPointId(2);
    cells->InsertNextCell(3, faceId);
    }
  delete [] faceId;

  polyData->SetPoints(this->GetAlignedPoints(index));
  polyData->SetPolys(cells);
  cells->Delete();
  return polyData;
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::RemoveShape(unsigned int id)
{
  itkDebugMacro(<< "Removing Shape " << id);
  unsigned int sz = m_Landmarks.size();
  if (id < sz)
    {
      m_Landmarks[id]->Delete();
      m_Landmarks.erase(m_Landmarks.begin() + id);
    }
  else
    std::cout << "Id out of Range " << id << std::endl;
  if(m_Landmarks.size() != (sz-1))
    std::cout << "Failed to delete ?" << std::endl;
  this->SetValid(false);
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::RemoveAllShapes()
{
  itkDebugMacro(<< "removing all shapes " <<  m_Landmarks.size());
  int sz = m_Landmarks.size();
  for (int i = sz-1; i >= 0; i--)
    {
      m_Landmarks[i]->Delete();
      m_Landmarks.erase(m_Landmarks.begin() + i);
    }
  if(m_Landmarks.size() != 0)
    std::cout << "Failed to delete all items ?" << std::endl;
  this->SetValid(false);
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::RemoveProcrustesAlignedPoints()
{
  int sz = m_ProcrustesAlignedPoints.size();
  for (int i = sz-1; i >= 0; i--)
    {
    m_ProcrustesAlignedPoints[i]->Delete();
    m_ProcrustesAlignedPoints.erase(m_ProcrustesAlignedPoints.begin() + i);
    }
  if(m_ProcrustesAlignedPoints.size() != 0)
    std::cout << "Failed to delete all items ?" << std::endl;
  this->SetValid(false);
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::UpdateProcrustesAlignment() throw(itk::ExceptionObject)
{
  vtkProcrustesAlignmentFilter *procrustesAlign = vtkProcrustesAlignmentFilter::New();

  this->RemoveProcrustesAlignedPoints();

  // Use rigid body alignment
  if(m_PoseType == 1)
    procrustesAlign->GetLandmarkTransform()->SetModeToRigidBody();
  if(m_PoseType == 2)
    procrustesAlign->GetLandmarkTransform()->SetModeToSimilarity();
  else if(m_PoseType == 3)
    procrustesAlign->GetLandmarkTransform()->SetModeToAffine();
  else
    procrustesAlign->GetLandmarkTransform()->SetModeToRigidBody();

  int sz = m_Landmarks.size();
#if VTK_MAJOR_VERSION <= 5
  procrustesAlign->SetNumberOfInputs(sz);
#endif
  procrustesAlign->StartFromCentroidOn();
  for(int i = 0; i < sz; i++)
    {
#if VTK_MAJOR_VERSION <= 5
    procrustesAlign->SetInput(i, this->GetShape(i));
#else
    procrustesAlign->AddInputDataObject(i, this->GetShape(i));
#endif
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
    vtkPolyData * aligned = vtkPolyData::New();
    vtkPoints * points = vtkPoints::New();
    points->DeepCopy(vtkPolyData::SafeDownCast(procrustesAlign->GetOutput(i))->GetPoints());
    aligned->SetPoints(points);
    points->Delete();
    m_ProcrustesAlignedPoints.push_back(aligned);
    }
  procrustesAlign->Delete();

  itkDebugMacro(<< "Performed Procrustes Alignment with " << sz << " shapes");
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurfaceSimilarityParameters(vtkPolyData * surface, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b, float bConstrain, unsigned int maxModesCount)
{
  this->GetSurfaceParametersToMatchShape1(surface, scale, tx, ty, tz, theta, phi, psi, b, bConstrain, maxModesCount);
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurfaceSimilarityParameters(vtkPolyData * surface, double &s, PointType &t, PointType &orientation, vtkFloatArray *b, float bConstrain, unsigned int maxModesCount)
{
  this->GetSurfaceSimilarityParameters(surface, s, t[0], t[1], t[2], orientation[0], orientation[1], orientation[2], b, bConstrain, maxModesCount);
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurfaceParametersToMatchShape1(vtkPolyData * surface, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b, float bConstrain, unsigned int maxModesCount)
{
  // Given a surface find the parameters that minimises the distance between the surface and the generated surfaceArea

  bool converged = false;
  int maximumIteration = 50;
  vtkPolyData * generatedSurface = NULL;
  vtkPolyData * oldGeneratedSurface = NULL;
  int i = 0;

  float * value = new float[1];
  for(int j = 0; j < b->GetNumberOfTuples(); j++)
    {
    value[0] = 0;
    b->SetTuple(j, value);
    }
  delete [] value;
  theta = 0; phi = 0; psi = 0;

  while((converged == false) && (i < maximumIteration))
    {
    // Given surface and b -> return pose
    this->GetSurfacePose(surface, scale, tx, ty, tz, theta, phi, psi, b);
    // Given surface, filter out pose -> return b
    this->GetSurfaceShapeParams(surface, scale, tx, ty, tz, theta, phi, psi, b);

    // restrict shape params to +- 3
    double * value2 = new double[1];
    for(int j = 0; j < b->GetNumberOfTuples(); j++)
      {
      b->GetTuple(j, value2);
      if(value2[0] > bConstrain)
        value2[0] = bConstrain;
      else if(value2[0] < -bConstrain)
        value2[0] = -bConstrain;
      b->SetTuple(j,value2);
      }
    delete [] value2;

    if(maxModesCount>0 && maxModesCount < b->GetNumberOfTuples())
    {
      value2[0] = 0.0;
      for(int j = maxModesCount; j < b->GetNumberOfTuples(); j++)
      {
        b->SetTuple(j,value2);
      }
    }

    generatedSurface = this->GetSurface(scale, tx, ty, tz, theta, phi, psi, b);

    if(i > 0)
      {
      converged = this->CheckConvergence(generatedSurface, oldGeneratedSurface, 1e-6);
      }
    if(oldGeneratedSurface != NULL)
      oldGeneratedSurface->Delete();
    oldGeneratedSurface = vtkPolyData::New();
    oldGeneratedSurface->DeepCopy(generatedSurface);
    generatedSurface->Delete();

    i++;
    }
  oldGeneratedSurface->Delete();
}

template <class TProfileSamplingPrecisionType>
bool
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::CheckConvergence(vtkPolyData * surface1, vtkPolyData * surface2, float minimumAverageSquareDistance)
{
  double score = this->PolyDataSquareDistance(surface1, surface2);
  int numberOfPoints = surface1->GetPoints()->GetNumberOfPoints();
  if(score/numberOfPoints < minimumAverageSquareDistance)
    return true;
  else
    return false;
}

template <class TProfileSamplingPrecisionType>
double
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::PolyDataDistance(vtkPolyData * surface1, vtkPolyData * surface2, vtkFloatArray* errorArray)
{
  double score = 0;
  int sz = surface1->GetPoints()->GetNumberOfPoints();
  vtkPoints * points1 = surface1->GetPoints();
  vtkPoints * points2 = surface2->GetPoints();

  if(errorArray == NULL)
    {
    for(int i = 0; i < sz; i++)
      {
      double * position1 = points1->GetPoint(i);
      double * position2 = points2->GetPoint(i);
      score += std::sqrt(vtkMath::Distance2BetweenPoints(position1, position2));
      }
    }
  else
    {
    for(int i = 0; i < sz; i++)
      {
      double * position1 = points1->GetPoint(i);
      double * position2 = points2->GetPoint(i);
      double value = std::sqrt(vtkMath::Distance2BetweenPoints(position1, position2));
      score += value;
      double * ScalarValue = new double[1];
      ScalarValue[0] = value;
      errorArray->SetTuple(i, ScalarValue);
      delete [] ScalarValue;
      }
    }

  return score;
}

template <class TProfileSamplingPrecisionType>
double
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::PolyDataSquareDistance(vtkPolyData * surface1, vtkPolyData * surface2, vtkFloatArray* errorArray)
{
  double score = 0;
  int sz = surface1->GetPoints()->GetNumberOfPoints();
  vtkPoints * points1 = surface1->GetPoints();
  vtkPoints * points2 = surface2->GetPoints();

  if(errorArray == NULL)
    {
    for(int i = 0; i < sz; i++)
      {
      double * position1 = points1->GetPoint(i);
      double * position2 = points2->GetPoint(i);
      score += vtkMath::Distance2BetweenPoints(position1, position2);
      }
    }
  else
    {
    for(int i = 0; i < sz; i++)
      {
      double * position1 = points1->GetPoint(i);
      double * position2 = points2->GetPoint(i);
      double value = vtkMath::Distance2BetweenPoints(position1, position2);
      score += value;
      double * ScalarValue = new double[1];
      ScalarValue[0] = value;
      errorArray->SetTuple(i, ScalarValue);
      delete [] ScalarValue;
      }
    }
  return score;
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
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
  if(b == NULL)
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

  for(int i = 0; i < 3; i++)
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
  for(int i = 0; i < 4; i++)
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
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurfacePose(vtkPolyData * surfaceX, double &scale, double &tx, double &ty, double &tz, double &theta, double &phi, double &psi, vtkFloatArray *b)
{
  // Translation is simply the centroid of the surface
  PointType centroidX = this->GetCentroid(surfaceX);
  tx = centroidX[0];
  ty = centroidX[1];
  tz = centroidX[2];
  //std::cout << "GetSurfacePose: Translation (" << tx << "," << ty << "," << tz << ")" << std::endl;

  // The optimal scaling factor is simply derived from the centroid size
  double scaleX = this->GetCentroidSize(surfaceX, centroidX);
  //std::cout << "GetSurfacePose: ScaleX " << scaleX << std::endl;

  // Generate shape we are aligning with
  vtkPolyData * surfaceY;
  if(b == NULL)
    {
    surfaceY = vtkPolyData::New();
    surfaceY->DeepCopy(this->GetMeanShape());
    }
  else
    {
    surfaceY = this->GetParameterisedShape(b);
    }
  PointType centroidY = this->GetCentroid(surfaceY);
  double scaleY = this->GetCentroidSize(surfaceY, centroidY);
  //cout << "GetSurfacePose: ScaleY " << scaleY << endl;

  scale = scaleX/scaleY;
  //cout << "GetSurfacePose: Scale " << scale << endl;

  // Filter out translation and scale effects
  vtkTransform * transformX = vtkTransform::New();
  transformX->PostMultiply();
  transformX->Identity();
  transformX->Translate(-centroidX[0], -centroidX[1], -centroidX[2]);
  transformX->Scale(1/scale, 1/scale, 1/scale);
  // Centroid Match
  vtkTransformPolyDataFilter *transformSurfaceX = vtkTransformPolyDataFilter::New();
  transformSurfaceX->SetTransform(transformX);
  transformSurfaceX->SetInputData(surfaceX);
  transformSurfaceX->Update();

  vtkTransform * transformY = vtkTransform::New();
  transformY->PostMultiply();
  transformY->Identity();
  transformY->Translate(-centroidY[0], -centroidY[1], -centroidY[2]);

  vtkTransformPolyDataFilter *transformSurfaceY = vtkTransformPolyDataFilter::New();
  transformSurfaceY->SetTransform(transformY);
  transformSurfaceY->SetInputData(surfaceY);
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
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurfaceShapeParams(vtkPolyData * surface, double scale, double tx, double ty, double tz, double theta, double phi, double psi, vtkFloatArray *b)
{
  // Remove pose effects -
  vtkTransform *transform = vtkTransform::New();
/*  transform->PostMultiply();
  transform->Identity();
  transform->Translate(-tx, -ty, -tz);
  transform->Scale(scale, scale, scale);
  transform->RotateZ(theta);
  transform->RotateX(phi);
  transform->RotateZ(psi);
  transform->Inverse();*/

  // Reverse transforms

  // Place
  vtkMatrix4x4 * matrix = this->GetSimilarityMatrix(scale, tx, ty, tz, theta, phi, psi, 0, 0, 0);
  matrix->Invert();
  transform->SetMatrix(matrix);
  matrix->Delete();
  //vtkMatrix4x4 * matrix = transform->GetMatrix();
  //cout << "Use Matrix " << endl;
  //cout << matrix->GetElement(0,0) << "," << matrix->GetElement(0,1) << "," << matrix->GetElement(0,2) << endl;
  //cout << matrix->GetElement(1,0) << "," << matrix->GetElement(1,1) << "," << matrix->GetElement(1,2) << endl;
  //cout << matrix->GetElement(2,0) << "," << matrix->GetElement(2,1) << "," << matrix->GetElement(2,2) << endl;

  vtkTransformPolyDataFilter *transformSurface = vtkTransformPolyDataFilter::New();
  transformSurface->SetTransform(transform);
  transformSurface->SetInputData(surface);
  transformSurface->Update();

  // Transformed surface should have zero origin,
  /*
  PointType centroid = this->GetCentroid(transformSurface->GetOutput());
  cout << "GetSurfaceShapeParams: Translation (" << centroid[0] << "," << centroid[1] << "," << centroid[2] << ")" << endl;
  cout << "GetSurfaceShapeParams: Scale " << this->GetCentroidSize(transformSurface->GetOutput()) << endl;

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(transformSurface->GetOutput());
  writer->SetFileName("GetSurfaceShapeParams.vtk");
  writer->Write();

  writer->SetInputData(surface);
  writer->SetFileName("GetSurfaceShapeParams_surface.vtk");
  writer->Write();
  writer->Delete();*/

  // Now find shape parameters
  vtkFloatArray * bNew = this->GetShapeParameters(transformSurface->GetOutput(), b->GetNumberOfTuples());

  float * value = new float[1];
  for(int i = 0; i < b->GetNumberOfTuples(); i++)
    {
    bNew->GetTypedTuple(i, value);
    b->SetTuple(i, value);
    }
  // Valgrind claims this isn't cleaned up
  delete [] value;

  bNew->Delete();
  transform->Delete();
  transformSurface->Delete();
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurfacePoseFromTransformMatrix(vnl_matrix<double> matrix, PointType &orientation)
{
  vtkMatrix4x4 * transformMatrix = vtkMatrix4x4::New();

  int rows = matrix.rows();
  int columns = matrix.columns();
  for(int i = 0; i < rows; i++)
    {
    for(int j = 0; j < columns; j++)
      {
      transformMatrix->SetElement(i,j, matrix(i,j));
      }
    }

  this->GetSurfacePoseFromTransformMatrix(transformMatrix, orientation);
  transformMatrix->Delete();
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurfacePoseFromTransformMatrix(vtkMatrix4x4 * matrix, PointType &orientation)
{
  // Calculate theta
  double theta = 0;

  double phi = 0;
  if(matrix->GetElement(2,2) > 0.99999999)
    {
    phi = 0;
    }
  else if(matrix->GetElement(2,2) < -0.99999999)
    {
    phi = vtkMath::Pi();
    }
  else
    {
    phi = acos(matrix->GetElement(2,2));
    }
  double psi = 0;

  if((fabs(phi) == 0) || (phi == vtkMath::Pi()))
    {
    // Need to solve other equations ;(...
    // We let theta = 0 -> sin(theta + psi)/cos(theta + psi) = m(0,1)/m(1,1)
    // -> theta + phi = atan2(m(0,1), m(1,1))
    theta = atan2(cos(phi)*matrix->GetElement(1,0), cos(phi)*matrix->GetElement(1,1));
    if(theta < 0)
      {
      theta = 2*vtkMath::Pi() + theta;
      }
    }
  else
    {
    double sinPhi = sin(phi);
    theta = atan2(matrix->GetElement(2,0)/sinPhi, matrix->GetElement(2,1)/sinPhi);
    if(theta < 0)
      {
      theta = 2*vtkMath::Pi() + theta;
      }
    psi = atan2(matrix->GetElement(0,2)/sinPhi, -matrix->GetElement(1,2)/sinPhi);
    if(psi < 0)
      {
      psi = 2*vtkMath::Pi() + psi;
      }
    }

  orientation[0] = 180.0*theta/vtkMath::Pi();
  orientation[1] = 180.0*phi/vtkMath::Pi();
  orientation[2] = 180.0*psi/vtkMath::Pi();
}

template <class TProfileSamplingPrecisionType>
double
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetCentroidSize(vtkPolyData * surface)
{
  PointType centroid = this->GetCentroid(surface);
  return this->GetCentroidSize(surface, centroid);
}

template <class TProfileSamplingPrecisionType>
double
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetCentroidSize(vtkPolyData * surface, PointType centroid)
{
  // TODO: Using true centroid size... but doesn't match definition of scale from Horn so using here.
  PointType pt;
  pt[0] = 0;
  pt[1] = 0;
  pt[2] = 0;
  double centroidSize = 0;
  for (int i = 0; i < surface->GetPoints()->GetNumberOfPoints(); i++)
    {
    double *p = surface->GetPoints()->GetPoint(i);
    pt[0] = p[0] - centroid[0];
    pt[1] = p[1] - centroid[1];
    pt[2] = p[2] - centroid[2];
    centroidSize += pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2];
    }
  return sqrt(centroidSize);
}

template <class TProfileSamplingPrecisionType>
typename StatisticalShapeModelBase<TProfileSamplingPrecisionType>::PointType
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetCentroid(vtkPolyData * surface)
{
  // Centroid
  PointType centroid;
  centroid[0] = 0;
  centroid[1] = 0;
  centroid[2] = 0;
  int numberPoints = surface->GetPoints()->GetNumberOfPoints();

  // Calculate center of shape
  double point[3];
  for (int i = 0; i < numberPoints; i++)
    {
    surface->GetPoints()->GetPoint(i, point);
    centroid[0] += point[0];
    centroid[1] += point[1];
    centroid[2] += point[2];
    }
  centroid[0] /= (float)numberPoints;
  centroid[1] /= (float)numberPoints;
  centroid[2] /= (float)numberPoints;
  return centroid;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurface(vtkMatrix4x4 * matrixTransform, vtkFloatArray *b, bool centroid)
{
  vtkPolyData * poly = vtkPolyData::New();

  vtkPolyData * surface;
  if(b == NULL)
    surface = this->GetMeanShape();
  else
    surface = this->GetParameterisedShape(b);

  vtkTransform * transform = vtkTransform::New();
  if(centroid == true)
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
  transformSurface->SetInputData(surface);
  transformSurface->Update();

  poly->DeepCopy(transformSurface->GetOutput());
  transformSurface->Delete();
  transform->Delete();
  surface->Delete();

  return poly;
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurface(vnl_matrix<double> pose, vtkFloatArray *b, bool centroid)
{
  vtkMatrix4x4 * matrixTransform = vtkMatrix4x4::New();
  for(int i = 0; i < 4; i++)
    {
    for(int j = 0; j < 4; j++)
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
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
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
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSurface(double s, PointType t, PointType orientation, vtkFloatArray *b)
{
  vtkMatrix4x4 * matrix = this->GetSimilarityMatrix(s, t[0], t[1], t[2], orientation[0], orientation[1], orientation[2]);
  vtkPolyData * surface =  this->GetSurface(matrix, b, true);
  matrix->Delete();
  return surface;
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::ClipSurfaceBasedOnScalars(vtkPolyData *surface, vtkPolyData *clippedSurface, vtkIdTypeArray *ids, const float thresholdValue)
{
  //~ std::cerr << "Clipping Surface" << std::endl;
  const int numberOfPoints = surface->GetNumberOfPoints();

  ids->SetName("PointIDs");
  ids->SetNumberOfTuples(numberOfPoints);
  ids->SetNumberOfComponents(1);
  ids->FillComponent(0, 0.0);

  if(!surface->GetPointData()->GetScalars())
  {
//    itkExceptionMacro(<< "The Surface doesn't have scalars set. Clipping surface will not work. Exiting.");
    std::cerr << "The Surface doesn't have scalars set. Clipping surface will not work. Exiting." << std::endl;
    return exit(EXIT_FAILURE);
  }

//  std::cout << "Number of Points " << numberOfPoints << std::endl;
  //extract regions of the mesh by thresholding
  vtkSmartPointer<vtkThreshold> thresholdSurface = vtkSmartPointer<vtkThreshold>::New();
#if VTK_MAJOR_VERSION <= 5
  thresholdSurface->SetInput(surface);
#else
  thresholdSurface->SetInputData(surface);
#endif
  thresholdSurface->ThresholdByUpper(thresholdValue);
  //Mark member marked outside parts as 0.0

  //connect raw threshold output
  vtkSmartPointer<vtkGeometryFilter> region = vtkSmartPointer<vtkGeometryFilter>::New();
  region->SetInputConnection(thresholdSurface->GetOutputPort());
  region->Update();
  vtkSmartPointer<vtkPolyData> mesh = region->GetOutput();

  //save id info for building cov matrix later by RSSM
  //numberOfPoints != mesh->GetNumberOfPoints() here
  for(int j = 0; j < mesh->GetNumberOfPoints(); j ++)
  {
    vtkIdType index = surface->FindPoint(mesh->GetPoint(j));
    ids->SetTuple1(index, j);
  }
  mesh->GetPointData()->AddArray(ids);

  clippedSurface->DeepCopy(mesh);
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::RestoreClippedSurfaceCorrespondence(vtkPolyData *clippedSurface, vtkIdTypeArray *ids, vtkPolyData *clippedCorrespondSurface)
{
  //~ std::cerr << "Clipping with Correspondence" << std::endl;
  vtkPoints *newPoints = vtkPoints::New();
  newPoints->SetNumberOfPoints(ids->GetNumberOfTuples());
  for(int k = 0; k < ids->GetNumberOfTuples(); k ++)
  {
    vtkIdType index = ids->GetValue(k);
    //~ cerr << k << " ";
    newPoints->SetPoint(k, clippedSurface->GetPoint(index));
  }
  clippedCorrespondSurface->SetPoints(newPoints);
  newPoints->Delete();
}

template <class TProfileSamplingPrecisionType>
vtkMatrix4x4 *
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSimilarityMatrix(double s, double tx, double ty, double tz, double theta, double phi, double psi, double centroidX, double centroidY, double centroidZ)
{
  vtkTransform *transform = vtkTransform::New();
  transform->PostMultiply();
  transform->Identity();
  transform->Translate(-centroidX, -centroidY, -centroidZ);
  transform->Scale(s, s, s);
  transform->RotateZ(theta);
  transform->RotateX(phi);
  transform->RotateZ(psi);
  transform->Translate(centroidX+tx, centroidY+ty, centroidZ+tz);

  vtkMatrix4x4 * matrix = vtkMatrix4x4::New();
  matrix->DeepCopy(transform->GetMatrix());
  transform->Delete();
  return matrix;
}

template <class TProfileSamplingPrecisionType>
vtkMatrix4x4 *
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSimilarityMatrix(double s, PointType t, PointType orientation, PointType centroid)
{
  vtkMatrix4x4 * transformMatrix = this->GetSimilarityMatrix(s, t[0], t[1], t[2], orientation[0], orientation[1], orientation[2], centroid[0], centroid[1], centroid[2]);
  return transformMatrix;
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModelBase<TProfileSamplingPrecisionType>
::GetSimilarityMatrix(vnl_matrix<double> &matrix, double s, PointType t, PointType orientation, PointType centroid)
{
  vtkMatrix4x4 * transformMatrix = this->GetSimilarityMatrix(s, t[0], t[1], t[2], orientation[0], orientation[1], orientation[2], centroid[0], centroid[1], centroid[2]);

  matrix.set_size(4,4);
  matrix.fill(0);

  for(int i = 0; i < 4; i++)
    {
    for(int j = 0; j < 4; j++)
      {
      matrix(i,j) = transformMatrix->GetElement(i,j);
      }
    }

  transformMatrix->Delete();
}

}

#endif
