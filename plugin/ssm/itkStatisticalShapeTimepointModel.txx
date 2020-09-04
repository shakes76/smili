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

#ifndef __itkStatisticalShapeTimepointModel_txx
#define __itkStatisticalShapeTimepointModel_txx

#include "itkStatisticalShapeTimepointModel.h"
//VTK
#include <vtkProcrustesAlignmentFilter.h>
#include <vtkMultiBlockDataSet.h>

namespace itk
{


template <class TProfileSamplingPrecisionType>
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
//~ ::StatisticalShapeTimepointModel() : StatisticalShapeModel<TProfileSamplingPrecisionType>()
::StatisticalShapeTimepointModel()
{
  this->m_PoseType = 2;
  if(Superclass::m_PCA != NULL)
    Superclass::m_PCA->Delete();
  Superclass::m_PCA = vtkPCAAnalysisTimepointFilter::New();
  Superclass::m_Precision = 0.9;
  Superclass::m_MeanScale = 0;
  Superclass::m_MeanShape = NULL;
  Superclass::m_Valid = false;
  m_TimepointMode = false;
  Reset();
}

template <class TProfileSamplingPrecisionType>
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::~StatisticalShapeTimepointModel()
{
//  vtkPCAAnalysisTimepointFilter *PCA = vtkPCAAnalysisTimepointFilter::SafeDownCast(Superclass::m_PCA);
//  PCA->Delete();
  itkDebugMacro(<<"StatisticalShapeTimepointModel Destructor");
}

template <class TProfileSamplingPrecisionType>
bool
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::AddAlignedShape(vtkPolyData *aligned)
{
  this->SetValid(false);
  // TODO: Add check if shape has same number of points and vertices
  // TODO: Maybe should also check it doesn't already exist !
  this->m_ProcrustesAlignedPoints.push_back(aligned);
  return true;
}


template <class TProfileSamplingPrecisionType>
bool
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::LoadModel(const char * filename) throw(itk::ExceptionObject)
{
  itkDebugMacro(<< "Loading FULL SSM with FileName " << filename);
  cout << "Loading FULL SSM with FileName " << filename << endl;

  typedef float headerType; //float for backwards compatibility to standard format
  typedef double writeType; 

  vtkPCAAnalysisTimepointFilter *PCA = vtkPCAAnalysisTimepointFilter::SafeDownCast(Superclass::m_PCA);

  std::ifstream fin(filename, std::ios::binary);
  if (!fin.is_open())
  {
    itkExceptionMacro(<< " SSM3D File " << filename << " does not exist ");
    return false;
  }

  int headerSize = 4;
  headerType *header = new headerType[headerSize];
  fin.read((char*)header, headerSize*sizeof(headerType));

  int totalShapes = header[0];
  int totalDimensions = header[1];
  int totalPoints = header[1]*header[2];
  int totalShapesPoints = totalShapes*totalPoints;
  int totalAlignedPoints = totalShapes*totalPoints;
  int totalMeanPoints = totalPoints;
  int totalCellvalues = header[3];
  int totalEigenvalues = totalShapes;
  int totalEigenvectors = totalShapesPoints;

  int dataSize = totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues
                 + totalEigenvalues + totalEigenvectors;
  writeType *data = new writeType[dataSize];

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
  if (!fin)
  {
    delete [] header;
    delete [] data;
    fin.close();
    itkDebugMacro(<< "Tried reading Full SSM data from file. Trying Standard SSM format instead for " << filename);
    cerr << "Tried reading Full SSM data from file. Trying Standard SSM format instead for " << filename << endl;
    return Superclass::LoadModel(filename);
  }
  m_Loaded = true;
  m_PreAligned = true;

  vtkIdType *faceId = new vtkIdType[3];
  int shapesCount = 0, alignedCount = 0;
  // For each shape get PointSet
  itkDebugMacro(<< "Loading Shapes");
#if VTK_MAJOR_VERSION <= 5
  PCA->SetNumberOfInputs(totalShapes);
#endif
  double *meanVector;
  for (int i = 0; i < totalShapes; i++)
  {
    vtkPolyData * shape = vtkPolyData::New();
    // Add points to shape i
    //std::cout << "About to add Points to polyData" << std::endl;
    vtkPoints *points=vtkPoints::New();
    for (int j = 0; j < header[2]; j++)
    {
      // Store each point in array before adding to pointset
      writeType *value = new writeType[totalDimensions];
      for (int k = 0; k < totalDimensions; k++)
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
      writeType *value = new writeType[totalDimensions];
      for (int k = 0; k < totalDimensions; k++)
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
      meanVector = new double[totalMeanPoints];
      for (int j = 0; j < totalMeanPoints; j++)
      {
        meanVector[j] = data[totalShapesPoints + totalAlignedPoints + j];
      }

      PCA->SetMeanShape(meanVector);
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

    polys->Delete();
    //std::cout << "Number of points in shape " << i << " is " << shape->GetNumberOfPoints() << std::endl;
    this->AddShape(shape);
    this->AddAlignedShape(aligned);
  #if VTK_MAJOR_VERSION <= 5
    PCA->SetInput(i, aligned);
  #else
    PCA->AddInputDataObject(i, aligned);
  #endif
  }

  //Get Eigenvalues and vectors.
  //Get Eigenvalues
  vtkFloatArray *eigenVals = vtkFloatArray::New();
  eigenVals->SetNumberOfTuples(totalEigenvalues);
  eigenVals->SetNumberOfComponents(1);
  for (int j = 0; j < eigenVals->GetNumberOfTuples(); j ++)
  {
    eigenVals->SetValue(j, data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + j]);
//    cout << eigenVals->GetValue(j) << ", ";
  }
//  cout << endl;
  PCA->SetEvals(eigenVals);
  eigenVals->Delete();
  itkDebugMacro(<< "Loaded Eigenvalues: " << PCA->GetEvals()->GetNumberOfTuples());
//  cerr << "Loaded Eigenvalues: " << PCA->GetEvals()->GetNumberOfTuples() << endl;

  //Get Eigenvectors
  double *matrix = new double[totalPoints*totalShapes];
  double **vectors = new double *[totalPoints];
  for(int i = 0; i < totalPoints; i++) {
    vectors[i] = &matrix[i*totalShapes];
  }
//  cerr << "Eigenvectors: " << endl;
  for (int j = 0; j < totalEigenvectors; j ++)
  {
    matrix[j] = data[totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues + totalEigenvalues + j];
//    cerr << matrix[j] << ", ";
  }
//  cerr << endl;
  PCA->SetEvectors(vectors);
  itkDebugMacro(<< "Loaded Eigenvectors: " << totalPoints << "x" << totalShapes);
//  cerr << "Loaded Eigenvectors: " << totalPoints << "x" << totalShapes << endl;

  //Set the weights in the PCA object
  this->UpdatePCA();

  for (int j = 0; j < totalShapes; j++)
  {
    for (int i = 0; i < header[2]; i++) {
      double x = vectors[i*3  ][j];
      double y = vectors[i*3+1][j];
      double z = vectors[i*3+2][j];

      vtkPolyData::SafeDownCast(PCA->GetOutput(j))->GetPoints()->SetPoint(i, x, y, z);
    }
  }

  if(Superclass::m_MeanShape != NULL)
    Superclass::m_MeanShape->Delete();
  Superclass::m_MeanShape = vtkPolyData::New();
  Superclass::m_MeanShape->DeepCopy(this->GetShape(0));
  vtkFloatArray* b = vtkFloatArray::New();
  PCA->GetParameterisedShape(b, Superclass::m_MeanShape);
  b->Delete();
  Superclass::m_MeanScale = (double)(this->GetCentroidSize(Superclass::m_MeanShape));

  this->SetValid(true);
  this->m_WeightMode = false;
  this->m_RobustMode = false;
  this->m_StandardMode = true;
  itkDebugMacro(<< "Updated Model");

  delete [] faceId;
  delete [] data;
  delete [] header;
  itkDebugMacro(<< "Finished Reading SSM Data File");
  return true;
}

template <class TProfileSamplingPrecisionType>
bool
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::SaveModel(const char * filename)
{
  itkDebugMacro(<< "Saving FULL SSM with filename " << filename);
  cout << "Saving FULL SSM with filename " << filename << endl;
  /// We want to save the current model to file
  /// File format is header, shape (point) data, aligned (point) data, connectiveness (triangle) data
  typedef float headerType; //float for backwards compatibility to standard format
  typedef double writeType; 

  vtkPCAAnalysisTimepointFilter *PCA = vtkPCAAnalysisTimepointFilter::SafeDownCast(Superclass::m_PCA);

  int sz = this->m_Landmarks.size();
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

    int totalShapes = header[0];
    int totalPoints = header[1]*header[2];
    int totalShapesPoints = header[0]*header[1]*header[2];
    int totalAlignedPoints = header[0]*header[1]*header[2];
    int totalMeanPoints = header[1]*header[2];
    int totalCellvalues = header[3];
    int totalEigenvalues = totalShapes;
    int totalEigenvectors = totalShapesPoints;

    itkDebugMacro(<< "Header is " << header[0] << " " << header[1] << " " << header[2] << " " << header[3]);
    int sizeData = totalShapesPoints + totalAlignedPoints + totalMeanPoints + totalCellvalues
                          + totalEigenvalues + totalEigenvectors;
    writeType *writer = new writeType [sizeData];

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
    for (int j = 0; j < totalMeanPoints; j++)
    {
      writer[count] = PCA->GetMeanShape()[j];
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

    ///Extract Eigenvalues
    //cerr<< "Extracting Eigenvalues.";
    vtkFloatArray *eigenVals = PCA->GetEvals();
    for (int j = 0; j < eigenVals->GetNumberOfTuples(); j++)
    {
      writer[count] = eigenVals->GetValue(j);
      count++;
    }

    ///Extract Eigenvectors
    //cerr<< "Extracting Eigenvectors.";
    double **eigenVecs = PCA->GetEvectors();
    for (int j = 0; j < totalPoints; j ++)
    {
      for(int k = 0; k < totalShapes; k ++)
      {
        writer[count] = eigenVecs[j][k];
        count++;
      }
    }

    // Write data to file
    itkDebugMacro(<< "Writing Data to File.");
    //cerr<< "Writing Data to File.";
    fout.write((char *)(header),headerSize*sizeof(headerType));
    fout.write((char *)(writer),sizeData*sizeof(writeType));
    fout.close();
    delete [] header;
    delete [] writer;
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
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
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
    }
  else
    std::cout << "Id out of Range " << id << std::endl;
  if(this->m_Landmarks.size() != (sz-1))
    std::cout << "Failed to delete ?" << std::endl;
  this->SetValid(false);
  this->m_Loaded = false;
  if(this->m_PCA)
    this->m_PCA->Delete();
  this->m_PCA = vtkPCAAnalysisTimepointFilter::New(); //reset
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::GetParameterisedShape(vtkFloatArray *b)
{
  this->Update();
  vtkPCAAnalysisTimepointFilter *PCA = vtkPCAAnalysisTimepointFilter::SafeDownCast(Superclass::m_PCA);
  //std::cout << "Number of tuples in b is " << b->GetNumberOfTuples() << std::endl;
  // Create a set of pts the same size as the number of landmarks
  vtkPolyData * shape = vtkPolyData::New();
  shape->DeepCopy(this->GetMeanShape());

  PCA->GetParameterisedShape(b, shape);

  return shape;
}

template <class TProfileSamplingPrecisionType>
vtkFloatArray *
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::GetShapeParameters(vtkPolyData *shape, int bsize)
{
  this->Update();
  vtkPCAAnalysisTimepointFilter *PCA = vtkPCAAnalysisTimepointFilter::SafeDownCast(Superclass::m_PCA);

  vtkFloatArray * b = vtkFloatArray::New();
  if(bsize == -1)
    bsize = PCA->GetModesRequiredFor(this->GetPrecision());
  b->SetNumberOfValues(bsize);
  PCA->GetShapeParameters(shape, b, bsize);
  return b;
}

template <class TProfileSamplingPrecisionType>
float
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
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
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::Update()
{
  if (this->GetValid() == false)
    this->GenerateData();
//  else if (m_Loaded)
//    cout << "Already up-to-date from load." << endl;
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::GenerateData()
{
//  cerr << "Timepoint SSM GenerateData" << endl;
  vtkPCAAnalysisTimepointFilter *PCA = vtkPCAAnalysisTimepointFilter::SafeDownCast(Superclass::m_PCA);
  PCA->SetTimepointMode(m_TimepointMode);

  this->UpdateProcrustesAlignment();
  this->UpdatePCA();

  this->RemoveProcrustesAlignedPoints();

  if(Superclass::m_MeanShape != NULL)
    Superclass::m_MeanShape->Delete();
  Superclass::m_MeanShape = vtkPolyData::New();
  Superclass::m_MeanShape->DeepCopy(this->GetShape(0));
  vtkFloatArray* b = vtkFloatArray::New();
  PCA->GetParameterisedShape(b, Superclass::m_MeanShape);
  b->Delete();
  Superclass::m_MeanScale = (double)(this->GetCentroidSize(Superclass::m_MeanShape));

  this->SetValid(true);
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::UpdateProcrustesAlignment() throw(itk::ExceptionObject)
{
  vtkProcrustesAlignmentFilter *procrustesAlign = vtkProcrustesAlignmentFilter::New();
  int sz = this->m_Landmarks.size();

  this->RemoveProcrustesAlignedPoints();

  if(!m_PreAligned)
  {
      // Use rigid body alignment
    if(this->m_PoseType == 1)
      procrustesAlign->GetLandmarkTransform()->SetModeToRigidBody();
    if(this->m_PoseType == 2)
      procrustesAlign->GetLandmarkTransform()->SetModeToSimilarity();
    else if(this->m_PoseType == 3)
      procrustesAlign->GetLandmarkTransform()->SetModeToAffine();
    else
      procrustesAlign->GetLandmarkTransform()->SetModeToRigidBody();

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
      this->m_ProcrustesAlignedPoints.push_back(aligned);
      }
    procrustesAlign->Delete();
    itkDebugMacro(<< "Performed Procrustes Alignment with " << sz << " shapes");
  }
  else
  {
    for (int i = 0; i < sz; i++)
    {
      vtkPolyData * aligned = vtkPolyData::New();
      aligned->DeepCopy(this->GetShape(i));
      this->m_ProcrustesAlignedPoints.push_back(aligned);
    }
  }
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::UpdatePCA() throw(itk::ExceptionObject)
{
  itkDebugMacro(<< "Update PCA");
  int sz = this->m_ProcrustesAlignedPoints.size();
  vtkPCAAnalysisTimepointFilter *PCA = vtkPCAAnalysisTimepointFilter::SafeDownCast(Superclass::m_PCA);

  if (m_Loaded)
    PCA->PreLoadedOn();
  else
  {
#if VTK_MAJOR_VERSION <= 5
    PCA->SetNumberOfInputs(sz);
#endif
  //    PCA->DebugOn();
    for (int i = 0; i < sz; i++)
    {
#if VTK_MAJOR_VERSION <= 5
      PCA->SetInput(i, this->m_ProcrustesAlignedPoints[i]);
#else
      PCA->AddInputDataObject(i, this->m_ProcrustesAlignedPoints[i]);
#endif
    }
  }

  try
  {
    PCA->Update();
  }
  catch ( itk::ExceptionObject & err )
  {
    itkExceptionMacro(<< "Failed PCA " << err);
  }

  this->m_Modes = PCA->GetModesRequiredFor(this->m_Precision);
  itkDebugMacro("Performed PCA Analysis with " << this->m_Modes << " significant modes");
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::Reset()
{
  m_Loaded = false;
  m_PreAligned = false;
}

template <class TProfileSamplingPrecisionType>
vtkPCAAnalysisTimepointFilter*
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::GetPCA()
{
  vtkPCAAnalysisTimepointFilter *PCA = vtkPCAAnalysisTimepointFilter::SafeDownCast(Superclass::m_PCA);

  return PCA;
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
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
#if VTK_MAJOR_VERSION <= 5
  transformSurface->SetInput(surface);
#else
  transformSurface->SetInputData(surface);
#endif
  transformSurface->Update();

  // Transformed surface should have zero origin,
  /*
  PointType centroid = this->GetCentroid(transformSurface->GetOutput());
  cout << "GetSurfaceShapeParams: Translation (" << centroid[0] << "," << centroid[1] << "," << centroid[2] << ")" << endl;
  cout << "GetSurfaceShapeParams: Scale " << this->GetCentroidSize(transformSurface->GetOutput()) << endl;

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(transformSurface->GetOutput());
#else
  writer->SetInputData(transformSurface->GetOutput());
#endif
  writer->SetFileName("GetSurfaceShapeParams.vtk");
  writer->Write();

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(surface);
#else
  writer->SetInputData(surface);
#endif
  writer->SetFileName("GetSurfaceShapeParams_surface.vtk");
  writer->Write();
  writer->Delete();*/

  // Now find shape parameters
  vtkFloatArray * bNew = GetShapeParameters(transformSurface->GetOutput(), b->GetNumberOfTuples());

  float * value = new float[1];
  for(int i = 0; i < b->GetNumberOfTuples(); i++)
    {
    bNew->GetTupleValue(i, value);
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
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
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
    surfaceY = GetParameterisedShape(b);
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
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
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
    surfaceY = GetParameterisedShape(b);
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
vtkPolyData *
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::GetSurface(vtkMatrix4x4 * matrixTransform, vtkFloatArray *b, bool centroid)
{
  vtkPolyData * poly = vtkPolyData::New();

  vtkPolyData * surface;
  if(b == NULL)
    surface = this->GetMeanShape();
  else
    surface = GetParameterisedShape(b);

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
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
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
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
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
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::GetSurface(double s, PointType t, PointType orientation, vtkFloatArray *b)
{
  vtkMatrix4x4 * matrix = this->GetSimilarityMatrix(s, t[0], t[1], t[2], orientation[0], orientation[1], orientation[2]);
  vtkPolyData * surface =  this->GetSurface(matrix, b, true);
  matrix->Delete();
  return surface;
}

template <class TProfileSamplingPrecisionType>
vtkMatrix4x4 *
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
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
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
::GetSimilarityMatrix(double s, PointType t, PointType orientation, PointType centroid)
{
  vtkMatrix4x4 * transformMatrix = this->GetSimilarityMatrix(s, t[0], t[1], t[2], orientation[0], orientation[1], orientation[2], centroid[0], centroid[1], centroid[2]);
  return transformMatrix;
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeTimepointModel<TProfileSamplingPrecisionType>
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
