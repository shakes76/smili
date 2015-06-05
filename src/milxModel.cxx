/*=========================================================================
  Program: milxSMILI
  Module: milxModel
  Author: Shekhar Chandra
  Language: C++
  Created: 18 Apr 2011 12:03:00 EST

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
#include "milxModel.h"

//VTK
#include <vtkAppendPolyData.h>
#include <vtkLandmarkTransform.h>
#include <vtkProcrustesAlignmentFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkGraphToPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkTable.h>
#include <vtkTriangle.h>
//VTK - Filters
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricDecimation.h>
#include <vtkQuadricClustering.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkLandmarkTransform.h>
#include <vtkCurvatures.h>
#include <vtkGradientFilter.h>
#include <vtkAssignAttribute.h>
#include <vtkSelectPolyData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkRotationFilter.h>
#include <vtkReflectionFilter.h>
#include <vtkThreshold.h> //clip
#include <vtkGeometryFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkDelaunay2D.h>
#include <vtkExtractEdges.h>
#include <vtkTubeFilter.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkIdFilter.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataPointSampler.h>
#include <vtkTensorGlyph.h>
#include <vtkHedgeHog.h>
#include <vtkPolyDataNormals.h>
#include <vtkFeatureEdges.h>
#include <vtkStripper.h>
#include <vtkConnectivityFilter.h>
#include <vtkSimpleElevationFilter.h>
#include <vtkElevationFilter.h>
#include <vtkKMeansStatistics.h>
#include <vtkQuantizePolyDataPoints.h>
//VTK 5.10 +
/*#include <vtkPolyDataToReebGraphFilter.h>
#include <vtkReebGraphSimplificationFilter.h>
#include <vtkReebGraphSurfaceSkeletonFilter.h>*/
//VTK - Sources
#include <vtkGlyphSource2D.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
//VTK Imaging
#include <vtkImageFlip.h>
#include <vtkImageMarchingCubes.h>
#if VTK_MAJOR_VERSION > 5
  #include <vtkMultiBlockDataSet.h>
  #include <vtkMultiBlockDataGroupFilter.h>
#endif
//VTK Extension
//#include "vtkAreaSimplificationMetric.h"
//SMILI
#include "milxMath.h"

typedef vtkDoubleArray milxArrayType;

namespace milx
{

Model::Model()
{
  InternalInPlaceOperation = false;
}

Model::Model(vtkSmartPointer<vtkPolyData> model)
{
  InternalInPlaceOperation = false;
  SetInput(model);
}

//IO
void Model::AddInput(vtkSmartPointer<vtkPolyData> model)
{
  vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
#if (VTK_MAJOR_VERSION > 5)
  if(CurrentModel)
    appendPolyData->AddInputData(CurrentModel);
    appendPolyData->AddInputData(model);
#else
  if(CurrentModel)
    appendPolyData->AddInput(CurrentModel);
    appendPolyData->AddInput(model);
#endif
  if(VTKProgressUpdates->IsBeingObserved())
    appendPolyData->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    appendPolyData->Update();

  UpdateModelState(appendPolyData);
}

void Model::SetInput(vtkSmartPointer<vtkPolyData> model)
{
  if(!CurrentModel)
    CurrentModel = vtkSmartPointer<vtkPolyData>::New();

  if(!InternalInPlaceOperation && CurrentModel != model) //!< Collection members can bypass copy here
    CurrentModel->DeepCopy(model);
  else
    CurrentModel = model;

  InputModel = CurrentModel;
}

void Model::SetTransform(vtkSmartPointer<vtkAbstractTransform> newTransform)
{
  if(!IsCurrentModel())
    return;

  if(!CurrentTransform)
  {
    CurrentTransform = vtkSmartPointer<vtkTransform>::New();
    CurrentTransform->Identity();
  }

  vtkSmartPointer<vtkTransform> linearTransform = vtkTransform::SafeDownCast(newTransform);
  if(linearTransform)
  {
    //CurrentTransform->Concatenate(linearTransform->GetMatrix());
    CurrentTransform->SetMatrix(linearTransform->GetMatrix());
  }

  PrintDebug("Transforming Surface");
  vtkSmartPointer<vtkTransformPolyDataFilter> transformedMesh = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    transformedMesh->SetInput(CurrentModel);
  #else
    transformedMesh->SetInputData(CurrentModel);
  #endif
    transformedMesh->SetTransform(newTransform);
  if(VTKProgressUpdates->IsBeingObserved())
    transformedMesh->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    transformedMesh->Update();

  UpdateModelState(transformedMesh);
}

void Model::SetGraph(vtkSmartPointer<vtkMutableUndirectedGraph> graph)
{
  if(!CurrentModel)
    CurrentModel = vtkSmartPointer<vtkPolyData>::New();

  vtkSmartPointer<vtkGraphToPolyData> graphPolyData = vtkSmartPointer<vtkGraphToPolyData>::New();
  #if VTK_MAJOR_VERSION <= 5
    graphPolyData->SetInput(graph);
  #else
    graphPolyData->SetInputData(graph);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    graphPolyData->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    graphPolyData->Update();

  UpdateModelState(graphPolyData);
}

void Model::SetPoints(vtkSmartPointer<vtkPoints> modelPoints)
{
  if(!CurrentModel)
    CurrentModel = vtkSmartPointer<vtkPolyData>::New();

  CurrentModel->SetPoints(modelPoints);
}

void Model::SetPolys(vtkSmartPointer<vtkCellArray> modelPolys)
{
  if(!IsCurrentModel())
    return;

  CurrentModel->SetPolys(modelPolys);
}

void Model::SetVectors(vtkSmartPointer<vtkDataArray> modelVectors)
{
  if(!IsCurrentModel())
    return;

  CurrentModel->GetPointData()->SetVectors(modelVectors);
}

void Model::SetScalars(vtkSmartPointer<vtkDataArray> modelScalars)
{
  if(!IsCurrentModel())
    return;

  CurrentModel->GetPointData()->SetScalars(modelScalars);
}

void Model::ResetScalars()
{
  FillScalars(0.0);
}

void Model::FillScalars(const coordinateType value)
{
  if(CurrentModel->GetPointData()->GetScalars())
    CurrentModel->GetPointData()->GetScalars()->FillComponent(0, value);
  else
  {
    vtkSmartPointer<milxArrayType> scalars = vtkSmartPointer<milxArrayType>::New();
      scalars->SetName("Reset Scalars");
      scalars->SetNumberOfTuples(CurrentModel->GetNumberOfPoints());
      scalars->SetNumberOfComponents(1);
      scalars->FillComponent(0, value);
      CurrentModel->GetPointData()->SetScalars(scalars);
  }
}

void Model::ClearScalars(vtkSmartPointer<vtkPolyData> &mesh)
{
  if(!mesh)
    return;

#if (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION < 8)
  mesh->GetPointData()->SetActiveScalars(NULL);
#else
  ///Strip existing scalars
  const int n = mesh->GetPointData()->GetNumberOfArrays();
  for(int j = 0; j < n; j ++)
  {
    //Ignore vector arrays etc.
    if(mesh->GetPointData()->GetArray(n-1-j)->GetNumberOfComponents() == 1)
      mesh->GetPointData()->RemoveArray(n-1-j);
  }
#endif
}

void Model::RemoveScalars()
{
  vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->DeepCopy(CurrentModel); //Deep copy is to allow undo

  ClearScalars(mesh);

  UpdateModelStateFromModel(mesh);
}

void Model::RemoveArrays()
{
  vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->DeepCopy(CurrentModel); //Deep copy is to allow undo

  ClearArrays(mesh);

  UpdateModelStateFromModel(mesh);
}

//Filters
void Model::Clean()
{
  if(!IsCurrentModel())
    return;

  PrintDebug("Cleaning");
  vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New(); ///Triangulate first
  #if VTK_MAJOR_VERSION <= 5
    cleanFilter->SetInput(CurrentModel);
  #else
    cleanFilter->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    cleanFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    cleanFilter->Update();

  UpdateModelState(cleanFilter);
}

void Model::Triangulate()
{
  if(!IsCurrentModel())
    return;

  PrintDebug("Triangulating");
  vtkSmartPointer<vtkTriangleFilter> triangulate = vtkSmartPointer<vtkTriangleFilter>::New(); ///Triangulate first
  #if VTK_MAJOR_VERSION <= 5
    triangulate->SetInput(CurrentModel);
  #else
    triangulate->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    triangulate->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    triangulate->Update();

  UpdateModelState(triangulate);
}

void Model::MatchInformation(vtkSmartPointer<vtkPolyData> matchMesh, const bool rescale)
{
  if(!IsCurrentModel())
    return;

  PrintDebug("Matching Information");
  coordinate modelCentroid = milx::Math<coordinateType>::Centroid(CurrentModel->GetPoints());
  coordinate otherCentroid = milx::Math<coordinateType>::Centroid(matchMesh->GetPoints());
#ifndef VTK_ONLY
  coordinate newCentroid = otherCentroid - modelCentroid;
#else
  coordinateType newCentroid[3];
  newCentroid[0] = otherCentroid[0] - modelCentroid[0];
  newCentroid[1] = otherCentroid[1] - modelCentroid[1];
  newCentroid[2] = otherCentroid[2] - modelCentroid[2];
#endif

  //Scale
  double modelScale = milx::Math<coordinateType>::CentroidSize(CurrentModel->GetPoints(), modelCentroid, true);
  double otherScale = milx::Math<coordinateType>::CentroidSize(matchMesh->GetPoints(), otherCentroid, true);
  double scaling = otherScale/modelScale;

  vtkSmartPointer<vtkTransform> transformer = vtkSmartPointer<vtkTransform>::New();
    transformer->Translate(newCentroid[0], newCentroid[1], newCentroid[2]);
    if(rescale)
        transformer->Scale(scaling, scaling, scaling);

  SetTransform(transformer);

#if defined VTK_ONLY
  delete modelCentroid;
  delete otherCentroid;
#endif
}

void Model::IterativeClosestPointsRegistration(vtkSmartPointer<vtkPolyData> fixedMesh, const bool similarity, const int maxPoints)
{
  if(!IsCurrentModel())
    return;

  PrintDebug("Iterative Closest Points");
  vtkSmartPointer<vtkIterativeClosestPointTransform> icp = vtkSmartPointer<vtkIterativeClosestPointTransform>::New();
    icp->SetSource(CurrentModel);
    icp->SetTarget(fixedMesh);
    icp->CheckMeanDistanceOn();
    icp->StartByMatchingCentroidsOn();
    icp->SetMaximumMeanDistance(0.001);
    icp->SetMaximumNumberOfIterations(300);
    icp->SetMaximumNumberOfLandmarks(maxPoints);
    icp->GetLandmarkTransform()->SetModeToRigidBody();
  if(VTKProgressUpdates->IsBeingObserved())
    icp->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    if(similarity)
      icp->GetLandmarkTransform()->SetModeToSimilarity();
    icp->Update();

  vtkMatrix4x4* transformMatrix = icp->GetMatrix();
  vtkSmartPointer<vtkTransform> transformer = vtkSmartPointer<vtkTransform>::New();
  transformer->SetMatrix(transformMatrix);

  SetTransform(transformer);
}

void Model::LandmarkBasedRegistration(vtkSmartPointer<vtkPolyData> fixedMesh, const bool similarity)
{
  if(!IsCurrentModel())
    return;

  if(CurrentModel->GetNumberOfPoints() != fixedMesh->GetNumberOfPoints())
  {
    PrintError("Landmark-based registration requires equal number of points. Ignoring Operation.");
    return;
  }

  PrintDebug("Landmark Transform");
  vtkSmartPointer<vtkLandmarkTransform> landmarkTransform = vtkSmartPointer<vtkLandmarkTransform>::New();
    landmarkTransform->SetModeToRigidBody();
    if(similarity)
      landmarkTransform->SetModeToSimilarity();
    landmarkTransform->SetSourceLandmarks(CurrentModel->GetPoints());
    landmarkTransform->SetTargetLandmarks(fixedMesh->GetPoints());
  if(VTKProgressUpdates->IsBeingObserved())
    landmarkTransform->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    landmarkTransform->Update();

  vtkMatrix4x4* transformMatrix = landmarkTransform->GetMatrix();
  vtkSmartPointer<vtkTransform> transformer = vtkSmartPointer<vtkTransform>::New();
  transformer->SetMatrix(transformMatrix);

  SetTransform(transformer);
}

void Model::Decimate(const float reduceFactor)
{
  Triangulate();
  ///Decimate using Quadric Decimation
  PrintDebug("PRO Decimation");
  vtkSmartPointer<vtkDecimatePro> decimateFilter = vtkSmartPointer<vtkDecimatePro>::New(); //Deliberate, destroy after use
  #if VTK_MAJOR_VERSION <= 5
    decimateFilter->SetInput(CurrentModel);
  #else
    decimateFilter->SetInputData(CurrentModel);
  #endif
    decimateFilter->SetTargetReduction(reduceFactor); //10% reduction (if there was 100 triangles, now there will be 90)
  if(VTKProgressUpdates->IsBeingObserved())
    decimateFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    decimateFilter->Update();

  UpdateModelState(decimateFilter);
}

void Model::QuadricDecimate(const float reduceFactor)
{
  Triangulate();
  ///Decimate using Quadric Decimation
  PrintDebug("Quadric Decimation");
  vtkSmartPointer<vtkQuadricDecimation> decimateFilter = vtkSmartPointer<vtkQuadricDecimation>::New(); //Deliberate, destroy after use
  #if VTK_MAJOR_VERSION <= 5
    decimateFilter->SetInput(CurrentModel);
  #else
    decimateFilter->SetInputData(CurrentModel);
  #endif
    decimateFilter->SetTargetReduction(reduceFactor);
  if(VTKProgressUpdates->IsBeingObserved())
    decimateFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    decimateFilter->Update();

  UpdateModelState(decimateFilter);
}

void Model::ClusterDecimate()
{
  Triangulate();
  ///Decimate using Quadric Decimation
  PrintDebug("Quadric Clustering Decimation");
  vtkSmartPointer<vtkQuadricClustering> decimateFilter = vtkSmartPointer<vtkQuadricClustering>::New(); //Deliberate, destroy after use
  #if VTK_MAJOR_VERSION <= 5
    decimateFilter->SetInput(CurrentModel);
  #else
    decimateFilter->SetInputData(CurrentModel);
  #endif
    decimateFilter->AutoAdjustNumberOfDivisionsOn();
  if(VTKProgressUpdates->IsBeingObserved())
    decimateFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    decimateFilter->Update();

  UpdateModelState(decimateFilter);
}

void Model::LaplacianSmoothing(const int iterations)
{
  if(!IsCurrentModel())
    return;

  PrintDebug("Laplacian Smoothing");
  vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New(); //Deliberate, destroy after use
  #if VTK_MAJOR_VERSION <= 5
    smoothFilter->SetInput(CurrentModel);
  #else
    smoothFilter->SetInputData(CurrentModel);
  #endif
    smoothFilter->SetNumberOfIterations(iterations);
  if(VTKProgressUpdates->IsBeingObserved())
    smoothFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    smoothFilter->Update();

  UpdateModelState(smoothFilter);
}

void Model::WindowedSincSmoothing(const int iterations, const float passBand)
{
  if(!IsCurrentModel())
    return;

  PrintDebug("Windowed Sinc Smoothing");
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoothFilter = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New(); //Deliberate, destroy after use
  #if VTK_MAJOR_VERSION <= 5
    smoothFilter->SetInput(CurrentModel);
  #else
    smoothFilter->SetInputData(CurrentModel);
  #endif
    smoothFilter->SetNumberOfIterations(iterations);
    smoothFilter->SetPassBand(passBand);
  if(VTKProgressUpdates->IsBeingObserved())
    smoothFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    smoothFilter->Update();

  UpdateModelState(smoothFilter);
}

void Model::Curvature(const bool meanCurvature)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkCurvatures> curvatureFilter = vtkSmartPointer<vtkCurvatures>::New(); //Deliberate, destroy after use
  #if VTK_MAJOR_VERSION <= 5
    curvatureFilter->SetInput(CurrentModel);
  #else
    curvatureFilter->SetInputData(CurrentModel);
  #endif
    curvatureFilter->SetCurvatureTypeToMean();
  if(!meanCurvature)
//      curvatureFilter->SetCurvatureTypeToGaussian();
      curvatureFilter->SetCurvatureTypeToMaximum();
//        curvatureFilter->SetCurvatureTypeToMinimum();
  if(VTKProgressUpdates->IsBeingObserved())
    curvatureFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    curvatureFilter->Update();

  UpdateModelState(curvatureFilter);
}

void Model::Gradient()
{
  if(!IsCurrentModel() || !CurrentModel->GetPointData()->GetScalars())
    return;

  vtkSmartPointer<vtkGradientFilter> gradientFilter = vtkSmartPointer<vtkGradientFilter>::New(); //Deliberate, destroy after use
  #if VTK_MAJOR_VERSION <= 5
    gradientFilter->SetInput(CurrentModel);
  #else
    gradientFilter->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    gradientFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive

  vtkSmartPointer<vtkGeometryFilter> geometry = vtkSmartPointer<vtkGeometryFilter>::New();
    geometry->SetInputConnection(gradientFilter->GetOutputPort());

  vtkSmartPointer<vtkAssignAttribute> vectors = vtkSmartPointer<vtkAssignAttribute>::New();
    vectors->SetInputConnection(geometry->GetOutputPort());
    vectors->Assign("Gradients", vtkDataSetAttributes::VECTORS, vtkAssignAttribute::POINT_DATA);
    vectors->Update();

  UpdateModelStateFromModel(vectors->GetPolyDataOutput());
}

void Model::DistanceMap(vtkSmartPointer<vtkPoints> loopPoints)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkSelectPolyData> selection = vtkSmartPointer<vtkSelectPolyData>::New();
  #if VTK_MAJOR_VERSION <= 5
    selection->SetInput(CurrentModel);
  #else
    selection->SetInputData(CurrentModel);
  #endif
    selection->SetLoop(loopPoints);
    selection->GenerateSelectionScalarsOn();
    selection->SetSelectionModeToSmallestRegion(); //negative scalars inside
//        selection->SetSelectionModeToClosestPointRegion();
  if(VTKProgressUpdates->IsBeingObserved())
    selection->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    selection->Update();

  UpdateModelState(selection);
}

void Model::GaussianWeightScalars(float stddev, float clampValue)
{
  if(!IsCurrentModel() || !GetScalars())
    return;

  vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->DeepCopy(CurrentModel);

  vtkSmartPointer<vtkDataArray> scalars = mesh->GetPointData()->GetScalars();

  const int n = mesh->GetNumberOfPoints();
  const float stddevSquared = 0.5*stddev*stddev;
  for(int j = 0; j < n; j ++)
  {
    double *value = scalars->GetTuple(j);
    value[0] = exp(-value[0]*value[0]/stddevSquared);
    if(value[0] < clampValue)
      value[0] = clampValue;
    scalars->SetTuple(j, value);
  }

  UpdateModelStateFromModel(mesh);
}

void Model::Rotate(const bool xAxis, const bool yAxis, const bool zAxis, const float angle, const coordinate centre)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkRotationFilter> rotator = vtkSmartPointer<vtkRotationFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    rotator->SetInput(CurrentModel);
  #else
    rotator->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    rotator->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    rotator->CopyInputOff();
    if(xAxis)
      rotator->SetAxisToX();
    if(yAxis)
      rotator->SetAxisToY();
    if(zAxis)
      rotator->SetAxisToZ();
    rotator->SetAngle(angle);
    rotator->SetCenter(centre[0], centre[1], centre[2]);
    rotator->SetNumberOfCopies(1);
    rotator->Update();

  ///connect raw threshold output
  vtkSmartPointer<vtkGeometryFilter> region = vtkSmartPointer<vtkGeometryFilter>::New();
    region->SetInputConnection(rotator->GetOutputPort());
  if(VTKProgressUpdates->IsBeingObserved())
    region->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    region->Update();

  UpdateModelState(region);
}

void Model::Flip(const bool xAxis, const bool yAxis, const bool zAxis)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkReflectionFilter> flipper = vtkSmartPointer<vtkReflectionFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    flipper->SetInput(CurrentModel);
  #else
    flipper->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    flipper->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    flipper->CopyInputOff();
    if(xAxis)
      flipper->SetPlaneToX();
    if(yAxis)
      flipper->SetPlaneToY();
    if(zAxis)
      flipper->SetPlaneToZ();
    flipper->Update();

  ///connect raw threshold output
  vtkSmartPointer<vtkGeometryFilter> region = vtkSmartPointer<vtkGeometryFilter>::New();
    region->SetInputConnection(flipper->GetOutputPort());
  if(VTKProgressUpdates->IsBeingObserved())
    region->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    region->Update();

  UpdateModelState(region);
}

void Model::Clip(const coordinateType belowValue, const coordinateType aboveValue)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkThreshold> thresh = vtkSmartPointer<vtkThreshold>::New();
  #if VTK_MAJOR_VERSION <= 5
    thresh->SetInput(CurrentModel);
  #else
    thresh->SetInputData(CurrentModel);
  #endif
    //~ thresh->ThresholdByLower(clusterID);
    //~ thresh->ThresholdByUpper(value);
    thresh->ThresholdBetween(belowValue, aboveValue);
  if(VTKProgressUpdates->IsBeingObserved())
    thresh->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive

  ///connect raw threshold output
  vtkSmartPointer<vtkGeometryFilter> region = vtkSmartPointer<vtkGeometryFilter>::New();
    region->SetInputConnection(thresh->GetOutputPort());
  if(VTKProgressUpdates->IsBeingObserved())
    region->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    region->Update();

  UpdateModelState(region);
}

void Model::Threshold(const coordinateType belowVal, const coordinateType aboveVal, const coordinateType outsideVal)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkPolyData> newMesh = vtkSmartPointer<vtkPolyData>::New();
    newMesh->DeepCopy(CurrentModel);

  //Threshold
  ThresholdScalars(newMesh, aboveVal, belowVal, outsideVal);

  UpdateModelStateFromModel(newMesh);
}

void Model::Threshold(const coordinateType belowVal, const coordinateType aboveVal, const coordinateType insideVal, const coordinateType outsideVal)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkPolyData> newMesh = vtkSmartPointer<vtkPolyData>::New();
    newMesh->DeepCopy(CurrentModel);

  //Threshold
  BinaryThresholdScalars(newMesh, aboveVal, belowVal, insideVal, outsideVal);

  UpdateModelStateFromModel(newMesh);
}

void Model::Mask(vtkSmartPointer<vtkPolyData> maskMesh)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkPolyData> newMesh = vtkSmartPointer<vtkPolyData>::New();
    newMesh->DeepCopy(CurrentModel);

  //Threshold
  MaskScalars(newMesh, maskMesh);

  UpdateModelStateFromModel(newMesh);
}

void Model::IsoSurface(vtkSmartPointer<vtkImageData> img, double minValue)
{
  vtkSmartPointer<vtkImageMarchingCubes> iso = vtkSmartPointer<vtkImageMarchingCubes>::New();
  #if VTK_MAJOR_VERSION <= 5
    iso->SetInput(img);
  #else
    iso->SetInputData(img);
  #endif
    iso->SetValue(0, minValue);
    iso->ComputeNormalsOn();
    iso->ComputeScalarsOff();
    //        iso->SetInputMemoryLimit(1000);
  if(VTKProgressUpdates->IsBeingObserved())
    iso->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    iso->Update();

  if(!IsCurrentModel())
    SetInput(iso->GetOutput());
  else
    UpdateModelState(iso);
}

vtkSmartPointer<vtkImageData> Model::Voxelise(const unsigned char insideValue, double *spacing, double *bounds, const size_t padVoxels)
{
  if(!IsCurrentModel())
    return NULL;

  vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();

  double totalBounds[6];
  if(bounds == NULL)
  {
      bounds = &totalBounds[0];
      CurrentModel->GetBounds(bounds);
  }
  whiteImage->SetSpacing(spacing);

  //Pad the bounds by a slice
  bounds[0] -= padVoxels*spacing[0];
  bounds[1] += padVoxels*spacing[0];
  bounds[2] -= padVoxels*spacing[1];
  bounds[3] += padVoxels*spacing[1];
  bounds[4] -= padVoxels*spacing[2];
  bounds[5] += padVoxels*spacing[2];

  // compute dimensions
  int dim[3];
  for (int i = 0; i < 3; i++)
  {
      dim[i] = static_cast<int>( ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]) );
  }
  PrintDebug("Dimensions of Voxelisation is " + NumberToString(dim[0]) + "x" + NumberToString(dim[1]) + "x" + NumberToString(dim[2]));
  whiteImage->SetDimensions(dim);
  whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

  double origin[3];
  // NOTE: I am not sure whether or not we had to add some offset!
  origin[0] = bounds[0];// + spacing[0] / 2;
  origin[1] = bounds[2];// + spacing[1] / 2;
  origin[2] = bounds[4];// + spacing[2] / 2;
  whiteImage->SetOrigin(origin);
#if VTK_MAJOR_VERSION <= 5
  whiteImage->SetNumberOfScalarComponents(1);
  whiteImage->SetScalarTypeToUnsignedChar();
  whiteImage->AllocateScalars();
#else
  whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif

  // fill the image with foreground voxels:
  unsigned char inval = insideValue;
  unsigned char outval = 0;
  for (vtkIdType i = 0; i < whiteImage->GetNumberOfPoints(); ++i)
      whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);

  // polygonal data --> image stencil:
  PrintDebug("Polygonal data --> Image stencil");
  vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
  #if VTK_MAJOR_VERSION <= 5
    pol2stenc->SetInput(CurrentModel);
  #else
    pol2stenc->SetInputData(CurrentModel);
  #endif
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
  if(VTKProgressUpdates->IsBeingObserved())
    pol2stenc->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    pol2stenc->Update();

  // cut the corresponding white image and set the background:
  PrintDebug("Image stencil --> Image");
  vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
  #if VTK_MAJOR_VERSION <= 5
    imgstenc->SetInput(whiteImage);
    imgstenc->SetStencil(pol2stenc->GetOutput());
  #else
    imgstenc->SetInputData(whiteImage);
    imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
  #endif
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(outval);
  if(VTKProgressUpdates->IsBeingObserved())
    imgstenc->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    imgstenc->Update();

  return imgstenc->GetOutput();
}

void Model::DelaunayGraph3D()
{
  if(!IsCurrentModel())
    return;

  ///Generate Delaunay Graph
  vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
  #if VTK_MAJOR_VERSION <= 5
    delaunay->SetInput(CurrentModel);
  #else
    delaunay->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    delaunay->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    delaunay->Update();
  /// Now we just pretty the mesh up with tubed edges and balls at the vertices.
  vtkSmartPointer<vtkExtractEdges> extract = vtkSmartPointer<vtkExtractEdges>::New();
    extract->SetInputConnection(delaunay->GetOutputPort());
  if(VTKProgressUpdates->IsBeingObserved())
    extract->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
  ///Tubes for the edges
  vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInputConnection(extract->GetOutputPort());
    tubes->SetRadius(0.02);
    tubes->SetNumberOfSides(3);
  if(VTKProgressUpdates->IsBeingObserved())
    tubes->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    tubes->Update();

  UpdateModelState(tubes);
}

void Model::Delaunay2DTriangulation()
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkCellArray> aCellArray = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyData> triangulatedModel = vtkSmartPointer<vtkPolyData>::New();
    triangulatedModel->SetPoints(CurrentModel->GetPoints());
    triangulatedModel->SetPolys(aCellArray);

  ///Generate Delaunay tessellation
  vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
  #if VTK_MAJOR_VERSION <= 5
    delaunay->SetInput(CurrentModel);
    delaunay->SetSource(triangulatedModel);
  #else
    delaunay->SetInputData(CurrentModel);
    delaunay->SetSourceData(triangulatedModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    delaunay->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    delaunay->Update();

  UpdateModelState(delaunay);
}

void Model::DelaunayTriangulation()
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
  #if VTK_MAJOR_VERSION <= 5
    delaunay->SetInput(CurrentModel);
  #else
    delaunay->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    delaunay->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    delaunay->Update();

  vtkSmartPointer<vtkGeometryFilter> region = vtkSmartPointer<vtkGeometryFilter>::New();
    region->SetInputConnection(delaunay->GetOutputPort());
  if(VTKProgressUpdates->IsBeingObserved())
    region->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    region->Update();

  UpdateModelState(region);
}

void Model::GenerateVertices()
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkVertexGlyphFilter> vertices = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    vertices->SetInput(CurrentModel);
  #else
    vertices->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    vertices->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    vertices->Update();

  UpdateModelState(vertices);
}

void Model::GenerateVertexScalars()
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkIdFilter> vertices = vtkSmartPointer<vtkIdFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    vertices->SetInput(CurrentModel);
  #else
    vertices->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    vertices->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    vertices->SetIdsArrayName("Point IDs");
    vertices->PointIdsOn();
    vertices->Update();

  UpdateModelStateFromModel(vertices->GetPolyDataOutput());
}

void Model::GenerateVerticesAs2DGlyphs(const GlyphType glyphType)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkGlyphSource2D> circle = vtkSmartPointer<vtkGlyphSource2D>::New();
  if(glyphType == Circle)
      circle->SetGlyphTypeToCircle();
  else if(glyphType == Triangle)
      circle->SetGlyphTypeToTriangle();
  else if(glyphType == Square)
      circle->SetGlyphTypeToSquare();
  else if(glyphType == Diamond)
      circle->SetGlyphTypeToDiamond();
  else
      circle->SetGlyphTypeToThickCross();
//           circle->SetScale(0.2);
//           circle->FilledOff();

  // Glyph the gradient vector
  vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
    glyph->SetSourceConnection(circle->GetOutputPort());
  #if VTK_MAJOR_VERSION <= 5
    glyph->SetInput(CurrentModel);
  #else
    glyph->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    glyph->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    glyph->Update();

  UpdateModelState(glyph);
}

void Model::GenerateTubes()
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    tubes->SetInput(CurrentModel);
  #else
    tubes->SetInputData(CurrentModel);
  #endif
    tubes->SetRadius(0.02);
    tubes->SetNumberOfSides(3);
  if(VTKProgressUpdates->IsBeingObserved())
    tubes->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    tubes->Update();

  UpdateModelState(tubes);
}

void Model::GeneratePointModel(const double newScale)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();

  vtkSmartPointer<vtkGlyph3D> glyphs = vtkSmartPointer<vtkGlyph3D>::New();
    glyphs->SetSourceConnection(sphere->GetOutputPort());
  #if VTK_MAJOR_VERSION <= 5
    glyphs->SetInput(CurrentModel);
  #else
    glyphs->SetInputData(CurrentModel);
  #endif
    glyphs->SetScaleModeToDataScalingOff();
    glyphs->SetScaleFactor(newScale);
  if(VTKProgressUpdates->IsBeingObserved())
    glyphs->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    glyphs->Update();

  UpdateModelState(glyphs);
}

void Model::GenerateSampledPoints(const double distance)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkPolyDataPointSampler> pointSampler = vtkSmartPointer<vtkPolyDataPointSampler>::New();
    pointSampler->SetDistance(distance);
  #if VTK_MAJOR_VERSION <= 5
    pointSampler->SetInput(CurrentModel);
  #else
    pointSampler->SetInputData(CurrentModel);
  #endif
    pointSampler->GenerateVertexPointsOn();
    pointSampler->GenerateVerticesOn();
  if(VTKProgressUpdates->IsBeingObserved())
    pointSampler->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    pointSampler->Update();

  UpdateModelState(pointSampler);
}

void Model::GenerateVectorField(const double newScale, const bool useNormals)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkArrowSource> arrows = vtkSmartPointer<vtkArrowSource>::New();
    arrows->SetShaftResolution(3);
    arrows->SetTipResolution(3);

  vtkSmartPointer<vtkGlyph3D> glyphs = vtkSmartPointer<vtkGlyph3D>::New();
    if(useNormals)
    {
      glyphs->SetVectorModeToUseNormal();
      glyphs->SetScaleModeToScaleByScalar();
    }
    else
    {
      glyphs->SetVectorModeToUseVector();
      glyphs->SetScaleModeToScaleByVector();
      glyphs->SetColorModeToColorByScalar();
    }
    glyphs->SetSourceConnection(arrows->GetOutputPort());
  #if VTK_MAJOR_VERSION <= 5
    glyphs->SetInput(CurrentModel);
  #else
    glyphs->SetInputData(CurrentModel);

  #endif
  if(newScale != 0.0)
    glyphs->SetScaleFactor(newScale);
  if(VTKProgressUpdates->IsBeingObserved())
    glyphs->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    glyphs->Update();

  UpdateModelState(glyphs);
}

void Model::GenerateTensorField(const double newScale)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetThetaResolution(4);
    sphere->SetPhiResolution(4);

  vtkSmartPointer<vtkTensorGlyph> ellipsoids = vtkSmartPointer<vtkTensorGlyph>::New();
  #if VTK_MAJOR_VERSION <= 5
    ellipsoids->SetInput(CurrentModel);
    ellipsoids->SetSource(sphere->GetOutput());
  #else
    ellipsoids->SetInputData(CurrentModel);
    ellipsoids->SetSourceData(sphere->GetOutput());
  #endif
    ellipsoids->ScalingOn();
    ellipsoids->ClampScalingOn();
    if(newScale != 0.0)
        ellipsoids->SetScaleFactor(newScale);
    ellipsoids->SetColorModeToScalars();
  if(VTKProgressUpdates->IsBeingObserved())
    ellipsoids->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    ellipsoids->Update();

  UpdateModelState(ellipsoids);
}

void Model::GenerateHedgehog(const double newScale)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkHedgeHog> hedgeHog = vtkSmartPointer<vtkHedgeHog>::New();
    hedgeHog->SetVectorModeToUseVector();
  #if VTK_MAJOR_VERSION <= 5
    hedgeHog->SetInput(CurrentModel);
  #else
    hedgeHog->SetInputData(CurrentModel);
  #endif
    if(newScale != 0.0)
      hedgeHog->SetScaleFactor(newScale);
  if(VTKProgressUpdates->IsBeingObserved())
    hedgeHog->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    hedgeHog->Update();

  UpdateModelState(hedgeHog);
}

void Model::GenerateNormals(const int pointNormals)
{
  if(!IsCurrentModel())
    return;

  PrintDebug("Generating Normals");
  vtkSmartPointer<vtkPolyDataNormals> surfaceNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
    surfaceNormals->ComputeCellNormalsOn();
    if(pointNormals == 0)
    {
      surfaceNormals->ComputePointNormalsOn();
      surfaceNormals->ComputeCellNormalsOff();
    }
    else if(pointNormals == 1)
    {
      surfaceNormals->ComputePointNormalsOff();
      surfaceNormals->ComputeCellNormalsOn();
    }
    else
    {
      surfaceNormals->ComputePointNormalsOn();
      surfaceNormals->ComputeCellNormalsOn();
    }
    surfaceNormals->SplittingOff();
    surfaceNormals->ConsistencyOn();
//    surfaceNormals->AutoOrientNormalsOn();
  #if VTK_MAJOR_VERSION <= 5
    surfaceNormals->SetInput(CurrentModel);
  #else
    surfaceNormals->SetInputData(CurrentModel);
  #endif
    surfaceNormals->SetFeatureAngle(37); //For shading
  if(VTKProgressUpdates->IsBeingObserved())
    surfaceNormals->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    surfaceNormals->Update();

  CurrentModel->GetPointData()->SetNormals( surfaceNormals->GetOutput()->GetPointData()->GetNormals() );
}

void Model::GenerateKMeansClustering(int numberOfClusters)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkTable> inputData = vtkSmartPointer<vtkTable>::New();
  for ( int c = 0; c < 3; ++c )
  {
    std::stringstream colName;
    colName << "coord " << c;
    vtkSmartPointer<vtkDoubleArray> doubleArray = vtkSmartPointer<vtkDoubleArray>::New();
      doubleArray->SetNumberOfComponents(1);
      doubleArray->SetName( colName.str().c_str() );
      doubleArray->SetNumberOfTuples(CurrentModel->GetNumberOfPoints());

    for ( int r = 0; r < CurrentModel->GetNumberOfPoints(); ++ r )
    {
        double p[3];
        CurrentModel->GetPoint(r, p);

        doubleArray->SetValue( r, p[c] );
    }

    inputData->AddColumn( doubleArray );
  }

  vtkSmartPointer<vtkKMeansStatistics> kMeansStatistics = vtkSmartPointer<vtkKMeansStatistics>::New();
  #if VTK_MAJOR_VERSION <= 5
    kMeansStatistics->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, inputData );
  #else
    kMeansStatistics->SetInputData( vtkStatisticsAlgorithm::INPUT_DATA, inputData );
  #endif
    kMeansStatistics->SetColumnStatus( inputData->GetColumnName( 0 ) , 1 );
    kMeansStatistics->SetColumnStatus( inputData->GetColumnName( 1 ) , 1 );
    kMeansStatistics->SetColumnStatus( inputData->GetColumnName( 2 ) , 1 );
    //kMeansStatistics->SetColumnStatus( "Testing", 1 );
    kMeansStatistics->RequestSelectedColumns();
    kMeansStatistics->SetAssessOption( true );
    kMeansStatistics->SetDefaultNumberOfClusters( numberOfClusters );
  if(VTKProgressUpdates->IsBeingObserved())
    kMeansStatistics->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    kMeansStatistics->Update();

  //Form scalars and add to model
  const float scaling = 256.0/numberOfClusters;
  vtkSmartPointer<vtkFloatArray> clusterArray = vtkSmartPointer<vtkFloatArray>::New();
    clusterArray->SetNumberOfComponents(1);
    clusterArray->SetName( "Cluster ID" );
  for (unsigned int r = 0; r < kMeansStatistics->GetOutput()->GetNumberOfRows(); r++)
  {
      vtkVariant v = kMeansStatistics->GetOutput()->GetValue(r,kMeansStatistics->GetOutput()->GetNumberOfColumns() - 1);
//        std::cout << "Point " << r << " is in cluster " << scaling*(v.ToInt()+1) << std::endl;
      clusterArray->InsertNextValue(scaling*(v.ToInt()+1));
  }
  CurrentModel->GetPointData()->SetScalars(clusterArray);
}

void Model::GenerateQuantisedPoints(float quantiseFactor)
{
  if(!IsCurrentModel())
    return;

  vtkSmartPointer<vtkQuantizePolyDataPoints> quantizeFilter = vtkSmartPointer<vtkQuantizePolyDataPoints>::New();
  #if VTK_MAJOR_VERSION <= 5
    quantizeFilter->SetInput(CurrentModel);
  #else
    quantizeFilter->SetInputData(CurrentModel);
  #endif
    quantizeFilter->SetQFactor(quantiseFactor);
  if(VTKProgressUpdates->IsBeingObserved())
    quantizeFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    quantizeFilter->Update();

  UpdateModelState(quantizeFilter);
}

void Model::GenerateCappedBoundaries()
{
  if(!IsCurrentModel())
    return;

  // Now extract feature edges, select the boundary
  vtkSmartPointer<vtkFeatureEdges> boundaryEdges = vtkSmartPointer<vtkFeatureEdges>::New();
  #if VTK_MAJOR_VERSION <= 5
    boundaryEdges->SetInput(CurrentModel);
  #else
    boundaryEdges->SetInputData(CurrentModel);
  #endif
    boundaryEdges->BoundaryEdgesOn();
    boundaryEdges->FeatureEdgesOff();
    boundaryEdges->NonManifoldEdgesOff();
    boundaryEdges->ManifoldEdgesOff();
  if(VTKProgressUpdates->IsBeingObserved())
    boundaryEdges->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive

  //Create triangles
  vtkSmartPointer<vtkStripper> boundaryStrips = vtkSmartPointer<vtkStripper>::New();
    boundaryStrips->SetInputConnection(boundaryEdges->GetOutputPort());
  if(VTKProgressUpdates->IsBeingObserved())
    boundaryStrips->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    boundaryStrips->Update();

  // Change the polylines into polygons
  vtkSmartPointer<vtkPolyData> boundaryPoly = vtkSmartPointer<vtkPolyData>::New();
    boundaryPoly->SetPoints(boundaryStrips->GetOutput()->GetPoints());
    boundaryPoly->SetPolys(boundaryStrips->GetOutput()->GetLines());

  AddInput(boundaryPoly);
}

void Model::GenerateRegions()
{
  if(!IsCurrentModel())
    return;

  // Now extract feature edges, select the boundary
  vtkSmartPointer<vtkConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkConnectivityFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    connectivityFilter->SetInput(CurrentModel);
  #else
    connectivityFilter->SetInputData(CurrentModel);
  #endif
    connectivityFilter->SetExtractionModeToAllRegions();
    connectivityFilter->ColorRegionsOn();
  if(VTKProgressUpdates->IsBeingObserved())
    connectivityFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive

  vtkSmartPointer<vtkGeometryFilter> region = vtkSmartPointer<vtkGeometryFilter>::New();
    region->SetInputConnection(connectivityFilter->GetOutputPort());
  if(VTKProgressUpdates->IsBeingObserved())
    region->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    region->Update();

  UpdateModelState(region);
}

void Model::GenerateElevation(double x, double y, double z)
{
  if(!IsCurrentModel())
    return;

  // Now extract feature edges, select the boundary
  vtkSmartPointer<vtkSimpleElevationFilter> elevationFilter = vtkSmartPointer<vtkSimpleElevationFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    elevationFilter->SetInput(CurrentModel);
  #else
    elevationFilter->SetInputData(CurrentModel);
  #endif
    elevationFilter->SetVector(x, y, z);
  if(VTKProgressUpdates->IsBeingObserved())
    elevationFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    elevationFilter->Update();

  UpdateModelStateFromModel(elevationFilter->GetPolyDataOutput());
}

void Model::GenerateElevation()
{
  if(!IsCurrentModel())
    return;

  double bounds[6];
  CurrentModel->GetBounds(bounds);

  // Now extract feature edges, select the boundary
  vtkSmartPointer<vtkElevationFilter> elevationFilter = vtkSmartPointer<vtkElevationFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    elevationFilter->SetInput(CurrentModel);
  #else
    elevationFilter->SetInputData(CurrentModel);
  #endif
    elevationFilter->SetLowPoint(bounds[0], bounds[2], bounds[4]);
    elevationFilter->SetHighPoint(bounds[0], bounds[2], bounds[5]);
  if(VTKProgressUpdates->IsBeingObserved())
    elevationFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    elevationFilter->Update();

  UpdateModelStateFromModel(elevationFilter->GetPolyDataOutput());
}
/*
void Model::GenerateReebGraph()
{
  if(!IsCurrentModel())
    return;

  // Now extract feature edges, select the boundary
  vtkSmartPointer<vtkPolyDataToReebGraphFilter> ReebGraphFilter = vtkSmartPointer<vtkPolyDataToReebGraphFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    ReebGraphFilter->SetInput(CurrentModel);
  #else
    ReebGraphFilter->SetInputData(CurrentModel);
  #endif
  if(VTKProgressUpdates->IsBeingObserved())
    ReebGraphFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    ReebGraphFilter->Update();

  // we don't define any custom simplification metric and use
  // the default one (persistence).
  vtkSmartPointer<vtkReebGraphSimplificationFilter> ReebSimplification = vtkSmartPointer<vtkReebGraphSimplificationFilter>::New();
//    vtkSmartPointer<vtkAreaSimplificationMetric> metric = vtkSmartPointer<vtkAreaSimplificationMetric>::New(); //vtk-ext
//    metric->SetLowerBound(0);
    // determining the maximum area
    //~ double globalArea = 0;
    //~ for(int i = 0; i < CurrentModel->GetNumberOfCells(); i ++)
    //~ {
      //~ vtkTriangle *t = vtkTriangle::SafeDownCast(CurrentModel->GetCell(i));
      //~ globalArea += t->ComputeArea();
    //~ }
    //~ metric->SetUpperBound(globalArea);
//    ReebSimplification->SetSimplificationMetric(metric);
    ReebSimplification->SetInputConnection(ReebGraphFilter->GetOutputPort()); //do not change, interface requires it
    ReebSimplification->SetSimplificationThreshold(0.01);
    ReebSimplification->Update();

  vtkSmartPointer<vtkReebGraphSurfaceSkeletonFilter> skeletonFilter = vtkSmartPointer<vtkReebGraphSurfaceSkeletonFilter>::New();
  #if VTK_MAJOR_VERSION <= 5
    skeletonFilter->SetInput(0, CurrentModel);
  #else
    skeletonFilter->SetInputData(0, CurrentModel);
  #endif
    skeletonFilter->SetInputConnection(1, ReebSimplification->GetOutputPort()); //do not change, interface requires it
    skeletonFilter->SetNumberOfSamples(5);
  if(VTKProgressUpdates->IsBeingObserved())
    skeletonFilter->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates); //Keeps UI responsive
    skeletonFilter->Update();
    vtkTable *surfaceSkeleton = skeletonFilter->GetOutput();

  PrintInfo("Reeb Graph has " + NumberToString(surfaceSkeleton->GetNumberOfColumns()*surfaceSkeleton->GetNumberOfRows()) + "points");
  vtkSmartPointer<vtkPolyData> embeddedSkeleton = vtkSmartPointer<vtkPolyData>::New();
    embeddedSkeleton->Allocate();

  vtkSmartPointer<vtkPoints> skeletonSamples = vtkSmartPointer<vtkPoints>::New();
    skeletonSamples->SetNumberOfPoints(surfaceSkeleton->GetNumberOfColumns()*surfaceSkeleton->GetNumberOfRows());

  int sampleId = 0;
  for(int i = 0; i < surfaceSkeleton->GetNumberOfColumns(); i++)
  {
    vtkDoubleArray *arc = vtkDoubleArray::SafeDownCast(surfaceSkeleton->GetColumn(i));

    // now add the samples to the skeleton polyData
    int initialSampleId = sampleId;
    for(int j = 0; j < arc->GetNumberOfTuples(); j++)
    {
      double point[3];
      arc->GetTupleValue(j, point);
      skeletonSamples->SetPoint(sampleId, point);
      sampleId++;
    }
    for(int j = 1; j < arc->GetNumberOfTuples(); j++)
    {
      vtkIdType samplePair[2];
      samplePair[0] = j - 1 + initialSampleId;
      samplePair[1] = j + initialSampleId;
      embeddedSkeleton->InsertNextCell(VTK_LINE, 2, samplePair);
    }
  }
  embeddedSkeleton->SetPoints(skeletonSamples);

  UpdateModelStateFromModel(embeddedSkeleton);
}*/

//Atomic
void Model::Append(vtkSmartPointer<vtkPolyData> model)
{
  vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
  #if VTK_MAJOR_VERSION <= 5
    if(CurrentModel) //Add current model too
      appendPolyData->AddInput(CurrentModel);
    appendPolyData->AddInput(model);
  #else
    if(CurrentModel) //Add current model too
      appendPolyData->AddInputData(CurrentModel);
    appendPolyData->AddInputData(model);
  #endif
    appendPolyData->Update();
    PreviousModel = appendPolyData->GetInput();
    CurrentModel = appendPolyData->GetOutput();
}

void Model::Scale(vtkSmartPointer<vtkPolyData> model, const coordinateType scale, const bool useCentroid)
{
  vtkSmartPointer<vtkPoints> points = model->GetPoints();
  const size_t numberOfPoints = model->GetNumberOfPoints();
#ifndef VTK_ONLY
  coordinate centroid(0.0);
#else
  coordinate centroid = new coordinateType[3];
#endif

  if(useCentroid)
    centroid = milx::Math<double>::Centroid(points);

  for(size_t j = 0; j < numberOfPoints; j ++)
  {
#ifndef VTK_ONLY
    coordinate point(points->GetPoint(j));

    if(useCentroid)
      point -= centroid;

    point *= scale;

    if(useCentroid)
      point += centroid;

    points->SetPoint(j, point.data_block());
#else
    coordinateType point[3];
    points->GetPoint(j, point);

    if(useCentroid)
    {
      point[0] -= centroid[0];
      point[1] -= centroid[1];
      point[2] -= centroid[2];
    }

    point[0] *= scale;
    point[1] *= scale;
    point[2] *= scale;

    if(useCentroid)
    {
      point[0] += centroid[0];
      point[1] += centroid[1];
      point[2] += centroid[2];
    }

    points->SetPoint(j, point);
#endif
  }
}

void Model::ScalarDifference(vtkSmartPointer<vtkPolyData> model1, vtkSmartPointer<vtkPolyData> model2, vtkSmartPointer<vtkPolyData> result)
{
  const size_t numberOfPoints = model1->GetNumberOfPoints();

  for(size_t j = 0; j < numberOfPoints; j ++)
  {
    double scalar1 = model1->GetPointData()->GetScalars()->GetTuple1(j);
    double scalar2 = model2->GetPointData()->GetScalars()->GetTuple1(j);
    double resultant = result->GetPointData()->GetScalars()->GetTuple1(j);

    resultant += scalar1 - scalar2;

    result->GetPointData()->GetScalars()->SetTuple1(j, resultant);
  }
}

void Model::ThresholdScalars(vtkSmartPointer<vtkPolyData> model, const coordinateType aboveVal, const coordinateType belowVal, const coordinateType outsideVal)
{
  BinaryThresholdScalars(model, aboveVal, belowVal, outsideVal, outsideVal);
}

void Model::BinaryThresholdScalars(vtkSmartPointer<vtkPolyData> model, const coordinateType aboveVal, const coordinateType belowVal, const coordinateType insideVal, const coordinateType outsideVal)
{
  if(!model->GetPointData()->GetScalars())
    return;

  vtkSmartPointer<vtkDataArray> scalars = model->GetPointData()->GetScalars();
  for(vtkIdType j = 0; j < scalars->GetNumberOfTuples(); j ++)
  {
    coordinateType scalar = scalars->GetTuple1(j);
    coordinateType value = 0.0;

    if(insideVal == outsideVal || vtkMath::IsNan(insideVal))
      value = scalar;
    else
      value = insideVal;

    if(scalar < belowVal || scalar > aboveVal)
      scalars->SetTuple1(j, outsideVal);
    else
      scalars->SetTuple1(j, value);
  }
}

void Model::MaskScalars(vtkSmartPointer<vtkPolyData> model1, vtkSmartPointer<vtkPolyData> model2)
{
  if(!model1->GetPointData()->GetScalars() || !model1->GetPointData()->GetScalars())
    return;

  vtkSmartPointer<vtkDataArray> scalars1 = model1->GetPointData()->GetScalars();
  vtkSmartPointer<vtkDataArray> scalars2 = model2->GetPointData()->GetScalars();
  for(vtkIdType j = 0; j < scalars1->GetNumberOfTuples(); j ++)
  {
    coordinateType scalar1 = scalars1->GetTuple1(j);
    coordinateType scalar2 = scalars2->GetTuple1(j);
    coordinateType value = 0.0;

    if(scalar2 >= 1)
      value = scalar1;

    scalars1->SetTuple1(j, value);
  }
}

vtkSmartPointer<vtkMatrix4x4> Model::ProcrustesAlign(vtkSmartPointer<vtkPolyData> source, vtkSmartPointer<vtkPolyData> target, bool rigid)
{
  if(source->GetNumberOfPoints() != target->GetNumberOfPoints())
    return NULL;

  double outPoint[3];
  vtkSmartPointer<vtkLandmarkTransform> landmarkTransform = vtkSmartPointer<vtkLandmarkTransform>::New();
    landmarkTransform->SetSourceLandmarks(source->GetPoints());
    landmarkTransform->SetTargetLandmarks(target->GetPoints());
    landmarkTransform->SetModeToSimilarity();
    if(rigid)
      landmarkTransform->SetModeToRigidBody();
    landmarkTransform->Update();
  for(vtkIdType v = 0; v < source->GetNumberOfPoints(); v ++)
  {
    landmarkTransform->InternalTransformPoint(source->GetPoint(v), outPoint);
    source->GetPoints()->SetPoint(v, outPoint);
  }

  vtkSmartPointer<vtkMatrix4x4> currentTransformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
    currentTransformMatrix->DeepCopy(landmarkTransform->GetMatrix());

  return currentTransformMatrix;
}

//Batch
void Model::ConcatenateCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
{
  const size_t n = collection->GetNumberOfItems();

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
    Append(collection->GetNextItem()); //!< concatenate the model
}

void Model::ScaleCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType scale)
{
  const size_t n = collection->GetNumberOfItems();

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
    Scale(collection->GetNextItem(), scale); //!< scale the model inplace
}

void Model::SmoothCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const size_t iterations)
{
  const size_t n = collection->GetNumberOfItems();
  InternalInPlaceOperation = true;

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();
    SetInput(mesh);
    WindowedSincSmoothing(iterations); //!< scale the model inplace
    mesh->DeepCopy(Result());
  }

  InternalInPlaceOperation = false;
}

void Model::LaplacianCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const size_t iterations)
{
  const size_t n = collection->GetNumberOfItems();
  InternalInPlaceOperation = true;

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();
    SetInput(mesh);
    LaplacianSmoothing(iterations); //!< scale the model inplace
    mesh->DeepCopy(Result());
  }

  InternalInPlaceOperation = false;
}

void Model::FlipCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const bool xAxis, const bool yAxis, const bool zAxis)
{
  const size_t n = collection->GetNumberOfItems();
  InternalInPlaceOperation = true;

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();
    SetInput(mesh);
    Flip(xAxis, yAxis, zAxis); //!< scale the model inplace
    mesh->DeepCopy(Result());
  }

  InternalInPlaceOperation = false;
}

void Model::DecimateCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType factor)
{
  const size_t n = collection->GetNumberOfItems();
  InternalInPlaceOperation = true;

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();
    SetInput(mesh);
    QuadricDecimate(factor); //!< scale the model inplace
    mesh->DeepCopy(Result());
  }

  InternalInPlaceOperation = false;
}

void Model::SplitCollection(vtkSmartPointer<vtkPolyDataCollection> collection, vtkSmartPointer<vtkPolyDataCollection> components, std::vector< vtkSmartPointer<vtkPolyDataCollection> > &splitCollections)
{
  const size_t n = collection->GetNumberOfItems();
  const size_t m = components->GetNumberOfItems();

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyDataCollection> splitCollection = vtkSmartPointer<vtkPolyDataCollection>::New();
    vtkSmartPointer<vtkPolyData> hybridSurface = collection->GetNextItem();
    vtkSmartPointer<vtkFloatArray> hybridScalars = vtkFloatArray::SafeDownCast(hybridSurface->GetPointData()->GetScalars());
    size_t accumulatedPoints = 0;

    components->InitTraversal();
    for(size_t k = 0; k < m; k ++)
    {
      vtkSmartPointer<vtkPolyData> splitResult = vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
      vtkSmartPointer<vtkPolyData> component = components->GetNextItem();
      const size_t numberOfPoints = component->GetNumberOfPoints();

      if(hybridScalars)
      {
        scalars->SetNumberOfComponents(hybridScalars->GetNumberOfComponents());
        scalars->SetNumberOfTuples(numberOfPoints);
      }
      
      //Split
      splitResult->DeepCopy(component); //Fast way to copy cell info etc.
      for(size_t l = 0; l < numberOfPoints; l ++) //We just change the geometric data
      {
        splitResult->GetPoints()->SetPoint(l, hybridSurface->GetPoint(l + accumulatedPoints));
        if(hybridScalars)
          scalars->SetTuple(l, hybridScalars->GetTuple(l + accumulatedPoints));
      }
      
      if(hybridScalars)
        splitResult->GetPointData()->SetScalars(scalars);

      splitCollection->AddItem(splitResult);
      accumulatedPoints += numberOfPoints;
    }

    splitCollections.push_back(splitCollection);
  }
}

void Model::ScalarDifferenceCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
{
  const size_t n = collection->GetNumberOfItems();

  collection->InitTraversal();
  vtkSmartPointer<vtkPolyData> surface1 = collection->GetNextItem();

  SetInput(surface1);
  ResetScalars();

  for(size_t j = 1; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> surface2 = collection->GetNextItem();

    ScalarDifference(surface1, surface2, CurrentModel); //!< threshold the model inplace
  }
}

void Model::ScalarThresholdCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType aboveVal, const coordinateType belowVal)
{
  const size_t n = collection->GetNumberOfItems();

  collection->InitTraversal();
  for(size_t j = 1; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> surface = collection->GetNextItem();

    ThresholdScalars(surface, aboveVal, belowVal, 0.0); //!< threshold the model inplace
  }
}

void Model::ClipCollection(vtkSmartPointer<vtkPolyDataCollection> collection, const coordinateType aboveVal, const coordinateType belowVal)
{
  const size_t n = collection->GetNumberOfItems();

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();
    SetInput(mesh);
    Clip(aboveVal, belowVal); //!< clip the model inplace
    mesh->DeepCopy(Result());
  }

  InternalInPlaceOperation = false;
}

void Model::ScalarDifferenceCollection(vtkSmartPointer<vtkPolyDataCollection> collection1, vtkSmartPointer<vtkPolyDataCollection> collection2, vtkSmartPointer<vtkPolyDataCollection> &results)
{
  const size_t n = collection1->GetNumberOfItems();

  if(collection1->GetNumberOfItems() != collection2->GetNumberOfItems())
  {
    PrintError("Unequal collections. Stopping");
    return;
  }

  collection1->InitTraversal();
  collection2->InitTraversal();
  for(size_t j = 0; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> surface1 = collection1->GetNextItem();
    vtkSmartPointer<vtkPolyData> surface2 = collection2->GetNextItem();

    vtkPolyData *result = vtkPolyData::New();
      result->DeepCopy(surface1);
    vtkSmartPointer<milxArrayType> scalars = vtkSmartPointer<milxArrayType>::New();
      scalars->SetNumberOfTuples(result->GetNumberOfPoints());
      scalars->SetNumberOfComponents(1);
      scalars->FillComponent(0, 0.0);
    result->GetPointData()->SetScalars(scalars);

    ScalarDifference(surface1, surface2, result); //!< diff the model scalars

    results->AddItem(result);
  }
}

void Model::ScalarStatisticsCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
{
  const size_t s = collection->GetNumberOfItems();

  collection->InitTraversal();
  vtkSmartPointer<vtkPolyData> firstSurface = collection->GetNextItem();
  const size_t n = firstSurface->GetNumberOfPoints();

  milxArrayType *minScalars = milxArrayType::New();
    minScalars->SetName("Minimum per Vertex");
    minScalars->SetNumberOfTuples(n);
    minScalars->SetNumberOfComponents(1);
    minScalars->FillComponent(0, minScalars->GetDataTypeMax());
  milxArrayType *maxScalars = milxArrayType::New();
    maxScalars->SetName("Maximum per Vertex");
    maxScalars->SetNumberOfTuples(n);
    maxScalars->SetNumberOfComponents(1);
    maxScalars->FillComponent(0, maxScalars->GetDataTypeMin());
  milxArrayType *meanScalars = milxArrayType::New();
    meanScalars->SetName("Mean per Vertex");
    meanScalars->SetNumberOfTuples(n);
    meanScalars->SetNumberOfComponents(1);
    meanScalars->FillComponent(0, 0.0);
  milxArrayType *varScalars = milxArrayType::New();
    varScalars->SetName("Variance per Vertex");
    varScalars->SetNumberOfTuples(n);
    varScalars->SetNumberOfComponents(1);
    varScalars->FillComponent(0, 0.0);

  size_t count = 0;
  collection->InitTraversal();
  for(size_t j = 0; j < s; j ++)
  {
    vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();

    if(!mesh->GetPointData()->GetScalars())
      continue; //ignore meshes with no scalars

    //safe casting requires double array
    vtkSmartPointer<vtkDataArray> scalars = mesh->GetPointData()->GetScalars();

    for(size_t k = 0; k < n; k ++)
    {
      if(scalars->GetTuple1(k) < minScalars->GetTuple1(k)) //min
        minScalars->SetTuple1(k, scalars->GetTuple1(k));
      if(scalars->GetTuple1(k) > maxScalars->GetTuple1(k)) //max
        maxScalars->SetTuple1(k, scalars->GetTuple1(k));

      meanScalars->SetTuple1( k, meanScalars->GetTuple1(k) + scalars->GetTuple1(k) ); //mean
    }

    count ++;
  }
  PrintDebug("Found "+NumberToString(count)+" Models with Scalars");

  for(size_t k = 0; k < n; k ++)
    meanScalars->SetTuple1( k, meanScalars->GetTuple1(k)/count ); //mean

  collection->InitTraversal();
  for(size_t j = 0; j < s; j ++)
  {
    vtkSmartPointer<vtkPolyData> mesh = collection->GetNextItem();

    if(!mesh->GetPointData()->GetScalars())
      continue; //ignore meshes with no scalars

    //safe casting requires double array
    vtkSmartPointer<vtkDataArray> scalars = mesh->GetPointData()->GetScalars();

    for(size_t k = 0; k < n; k ++)
    {
      const coordinateType diffValue = scalars->GetTuple1(k) - meanScalars->GetTuple1(k);
      varScalars->SetTuple1( k, varScalars->GetTuple1(k) + (diffValue*diffValue)/count ); //variance
    }
  }

  //Create resulting mesh
  SetInput(firstSurface);
  RemoveScalars();
  CurrentModel->GetPointData()->AddArray(meanScalars);
  CurrentModel->GetPointData()->AddArray(minScalars);
  CurrentModel->GetPointData()->AddArray(maxScalars);
  CurrentModel->GetPointData()->AddArray(varScalars);
  CurrentModel->GetPointData()->SetActiveScalars("Mean per Vertex");
}

void Model::ScalarRemoveCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
{
  const size_t n = collection->GetNumberOfItems();

  collection->InitTraversal();
  for(size_t j = 0; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> surface = collection->GetNextItem();

    ClearScalars(surface); //!< inplace
  }
}

void Model::ScalarCopyCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
{
  const size_t n = collection->GetNumberOfItems();

  collection->InitTraversal();
  vtkSmartPointer<vtkPolyData> templateSurface = collection->GetNextItem();
  for(size_t j = 1; j < n; j ++)
  {
    vtkSmartPointer<vtkPolyData> surface = collection->GetNextItem();

    ClearScalars(surface); //!< inplace

    surface->GetPointData()->SetScalars(templateSurface->GetPointData()->GetScalars());
  }
}

double Model::MeanSquaredErrorCollection(vtkSmartPointer<vtkPolyDataCollection> collection)
{
  const size_t n = collection->GetNumberOfItems();
  size_t count = 0;
  double error = 0.0;

  collection->InitTraversal();
  while(count < n)
  {
    vtkSmartPointer<vtkPolyData> surface1 = collection->GetNextItem();
    count ++;

    if(count >= n)
      break;

    vtkSmartPointer<vtkPolyData> surface2 = collection->GetNextItem();
    count ++;

    error += milx::Math<coordinateType>::MeanSquaredError(surface1->GetPoints(), surface2->GetPoints());
  }

  return error;
}

vtkSmartPointer<vtkPolyDataCollection> Model::ProcrustesAlignCollection(vtkSmartPointer<vtkPolyDataCollection> collection, bool rigid)
{
  const int n = collection->GetNumberOfItems();

  collection->InitTraversal();

  //Alignment filter
  vtkSmartPointer<vtkProcrustesAlignmentFilter> procrustes = vtkSmartPointer<vtkProcrustesAlignmentFilter>::New();
  #if VTK_MAJOR_VERSION <=5
    procrustes->SetNumberOfInputs(n);
  #endif
    if(rigid)
    {
      procrustes->GetLandmarkTransform()->SetModeToRigidBody();
      procrustes->StartFromCentroidOn();
    }
    else
      procrustes->GetLandmarkTransform()->SetModeToSimilarity();

  //Current model used for mean

  //mean of scalars
  vtkSmartPointer<milxArrayType> meanScalars = vtkSmartPointer<milxArrayType>::New();
    meanScalars->SetNumberOfComponents(1);

#if VTK_MAJOR_VERSION <=5
  for(int j = 0; j < n; j ++)
    procrustes->SetInput(j, collection->GetNextItem());
#else
  vtkSmartPointer<vtkMultiBlockDataGroupFilter> group = vtkSmartPointer<vtkMultiBlockDataGroupFilter>::New();
  for(int j = 0; j < n; j ++)
    group->AddInputData(collection->GetNextItem());
  procrustes->SetInputConnection(group->GetOutputPort());
#endif

  int numberOfPoints = 0;
  if(n > 0)
  {
  #if VTK_MAJOR_VERSION <=5
    numberOfPoints = procrustes->GetInput(0)->GetNumberOfPoints();
  #else
    numberOfPoints = vtkPointSet::SafeDownCast(procrustes->GetOutput()->GetBlock(0))->GetNumberOfPoints();
  #endif
    meanScalars->SetNumberOfTuples(numberOfPoints);
    meanScalars->FillComponent(0, 0.0);
  }
  procrustes->Update();

  vtkSmartPointer<vtkPolyDataCollection> alignedCollection = vtkSmartPointer<vtkPolyDataCollection>::New();
  collection->InitTraversal();
  for(int j = 0; j < n; j ++)
  {
    vtkPolyData *aligned = vtkPolyData::New();
      aligned->DeepCopy(collection->GetNextItem());
    #if VTK_MAJOR_VERSION <=5
      aligned->SetPoints(procrustes->GetOutput(j)->GetPoints());
    #else
      aligned->SetPoints( vtkPointSet::SafeDownCast(procrustes->GetOutput()->GetBlock(j))->GetPoints() );
    #endif

    if(aligned->GetPointData()->GetScalars())
    {
      for(int k = 0; k < numberOfPoints; k ++)
      {
        const double scalar = aligned->GetPointData()->GetScalars()->GetTuple1(k);
        double resultant = meanScalars->GetTuple1(k);

        resultant += scalar;

        meanScalars->SetTuple1(k, resultant);
      }
    }

    if(j == 0)
      SetInput(aligned);

    alignedCollection->AddItem(aligned);
  }

  //Finalise mean mesh
  for(int k = 0; k < numberOfPoints; k ++)
  {
    double resultant = meanScalars->GetTuple1(k);

    resultant /= n; //norm

    meanScalars->SetTuple1(k, resultant);
  }
  CurrentModel->SetPoints(procrustes->GetMeanPoints());
  CurrentModel->GetPointData()->SetScalars(meanScalars);

  return alignedCollection;
}

vtkSmartPointer<vtkPolyDataCollection> Model::IterativeClosestPointsAlignCollection(vtkSmartPointer<vtkPolyDataCollection> collection, bool rigid,
    vtkSmartPointer<vtkTransformCollection> tCollection)
{
  const int n = collection->GetNumberOfItems();
  vtkSmartPointer<vtkPolyDataCollection> alignedCollection = vtkSmartPointer<vtkPolyDataCollection>::New();

  //Use the last item as the reference "fixed" surface
  vtkPolyData *last = static_cast<vtkPolyData *>(collection->GetItemAsObject(n-1));

  collection->InitTraversal();
  for(int j = 0; j < n-1; j ++)
  {
    vtkPolyData *aligned = collection->GetNextItem();

    SetInput(aligned);
    MatchInformation(last, !rigid);
    IterativeClosestPointsRegistration(last, !rigid);

    vtkSmartPointer<vtkPolyData> result = vtkSmartPointer<vtkPolyData>::New();
    result->DeepCopy(Result());
    alignedCollection->AddItem(result);

    if (tCollection.GetPointer() != 0) {
      vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
      transform->DeepCopy(this->CurrentTransform);
      tCollection->AddItem(transform);
      //tCollection->AddItem(vtkSmartPointer<vtkTransform>::NewInstance(this->CurrentTransform));
    }
  }

  return alignedCollection;
}

} //end milx namespace
