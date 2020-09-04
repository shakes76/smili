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

#ifndef __itkStatisticalShapeModel_txx
#define __itkStatisticalShapeModel_txx

#include "itkStatisticalShapeModel.h"

#include <fstream>
#include <iostream>

#include <vtkFloatArray.h>
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include <vtkMatrix4x4.h>

#include <vtkProcrustesAlignmentFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLandmarkTransform.h>
#include <vtkPolyDataWriter.h>

namespace itk
{

/**
 *
 */
template <class TProfileSamplingPrecisionType>
StatisticalShapeModel<TProfileSamplingPrecisionType>
::StatisticalShapeModel()
{
  m_PCA = vtkPCAAnalysisFilter::New();
}

template <class TProfileSamplingPrecisionType>
StatisticalShapeModel<TProfileSamplingPrecisionType>
::~StatisticalShapeModel()
{
  if(m_PCA != NULL)
    m_PCA->Delete();
  itkDebugMacro(<<"StatisticalShapeModel Destructor");
}

template <class TProfileSamplingPrecisionType>
void
StatisticalShapeModel<TProfileSamplingPrecisionType>
::GenerateData()
{
  this->UpdateProcrustesAlignment();
  this->UpdatePCA();

  this->RemoveProcrustesAlignedPoints();

  if(this->m_MeanShape != NULL)
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
vtkFloatArray *
StatisticalShapeModel<TProfileSamplingPrecisionType>
::GetShapeParameters(vtkPolyData *shape, int bsize)
{
  this->Update();

  vtkFloatArray * b = vtkFloatArray::New();
  if(bsize == -1)
    bsize = this->m_PCA->GetModesRequiredFor(this->GetPrecision());
  b->SetNumberOfValues(bsize);
  m_PCA->GetShapeParameters(shape, b, bsize);
  return b;
}

template <class TProfileSamplingPrecisionType>
float
StatisticalShapeModel<TProfileSamplingPrecisionType>
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
StatisticalShapeModel<TProfileSamplingPrecisionType>
::UpdatePCA() throw(itk::ExceptionObject)
{
  int sz = this->m_ProcrustesAlignedPoints.size();
#if VTK_MAJOR_VERSION <= 5
  m_PCA->SetNumberOfInputs(sz);
#endif
  for(int i = 0; i < sz; i++)
    {
#if VTK_MAJOR_VERSION <= 5
    m_PCA->SetInput(i, this->m_ProcrustesAlignedPoints[i]);
#else
    m_PCA->AddInputDataObject(i, this->m_ProcrustesAlignedPoints[i]);
#endif
    }
  try
    {
    m_PCA->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    itkExceptionMacro(<< "Failed PCA " << err);
    }
  this->m_Modes = m_PCA->GetModesRequiredFor(this->m_Precision);
  itkDebugMacro("Performed PCA Analysis with " << this->m_Modes << " significant modes");
}

template <class TProfileSamplingPrecisionType>
vtkPolyData *
StatisticalShapeModel<TProfileSamplingPrecisionType>
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
double
StatisticalShapeModel<TProfileSamplingPrecisionType>
::MahalanobisDistanceBetween(vtkPolyData * poly1, vtkPolyData * poly2, const bool includeScaling)
{
  //Setup eigenfeature vectors
  vtkFloatArray *b1 = vtkFloatArray::New();
    b1->SetNumberOfComponents(1);
    b1->SetNumberOfTuples(this->GetNumberOfModes());
  vtkFloatArray *b2 = vtkFloatArray::New();
    b2->SetNumberOfComponents(1);
    b2->SetNumberOfTuples(this->GetNumberOfModes());

  // Now find shape parameters by stripping pose parameters and projecting shapes into model
  double scale1, scale2, tx, ty, tz, theta, phi, psi; //Unused
  this->GetSurfaceSimilarityParameters(poly1, scale1, tx, ty, tz, theta, phi, psi, b1);
  this->GetSurfaceSimilarityParameters(poly2, scale2, tx, ty, tz, theta, phi, psi, b2);

  double overallScaling = 1.0;
  if(includeScaling)
  {
    if(scale2 < scale1)
      overallScaling = scale2/scale1;
    else
      overallScaling = scale1/scale2;
  }

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
    distance += overallScaling*b1->GetValue(i)*b2->GetValue(i)/sqrt(eigenvalues->GetValue(i));
//    distance += b1->GetValue(i)*b2->GetValue(i);
  distance /= norm1*norm2; //Norm

  b1->Delete();
  b2->Delete();

  return -distance;
}

template <class TProfileSamplingPrecisionType>
bool
StatisticalShapeModel<TProfileSamplingPrecisionType>
::LoadModel(const char * filename) throw(itk::ExceptionObject)
{
  itkDebugMacro(<< "Loading SSM with FileName " << filename);
  std::ifstream fin(filename, std::ios::binary);
  if(!fin.is_open())
    {
    itkExceptionMacro(<< "SSM3D File " << filename << " does not exist ");
    return false;
    }

  int headerSize = 4;
  float *header = new float[headerSize];
  fin.read((char*)header, headerSize*sizeof(float));

  int dataSize = (int)(header[0]*header[1]*header[2] + header[3]);

  //std::cout << "Amount of data to read is " << dataSize << std::endl;
  itkDebugMacro(<< "Number of shapes is " << header[0]);
  itkDebugMacro(<< "Number of Dimensions " << header[1]);
  itkDebugMacro(<< "Number of Landmarks " << header[2]);
  itkDebugMacro(<< "Number of Faces " << header[3]/3.0);

  if((header[0] <=0) || (header[1] <= 0) || (header[2] <= 0) || (header[3] <=0))
    {
    delete [] header;
    fin.close();
    itkExceptionMacro(<< "Not a valid ssm " << filename);
    return false;
    }
  this->RemoveAllShapes();

  // Check file size to ensure it's valid using header info
  long begin,end;
  begin = fin.tellg();
  fin.seekg (0, std::ios::end);
  end = fin.tellg();

  if( (int)((end-begin)/sizeof(float)) != dataSize)
    {
    delete [] header;
    fin.close();
    itkExceptionMacro(<< "SSM file has length: " << (end-begin)/sizeof(float) << " which does not match expected dataSize=" << dataSize);
    return false;
    }

  fin.seekg(begin);
  float *data = new float[dataSize];
  fin.read((char*)data, dataSize*sizeof(float));
  if(!fin)
    {
    delete [] header;
    delete [] data;
    fin.close();
    itkExceptionMacro(<< "Failed reading SSM data from file " << filename);
    return false;
    }

  int count = 0;
  vtkIdType *faceId = new vtkIdType[3];
  // For each shape get PointSet
  for (int i = 0; i < header[0]; i++)
    {
    vtkPolyData * shape = vtkPolyData::New();
    // Add points to shape i
    //std::cout << "About to add Points to polyData" << std::endl;
    vtkPoints *points=vtkPoints::New();
    for (int j = 0; j < header[2]; j++)
      {
      // Store each point in array before adding to pointset
      float *value = new float[(int)(header[1])];
      for (int k = 0; k < header[1]; k++)
        {
        value[k] = data[count];
        count++;
        }
      points->InsertNextPoint(value);
      delete [] value;
      }
    //std::cout << "About to add Faces to polyData" << std::endl;
    // Add faces for shape i
    vtkCellArray *polys = vtkCellArray::New();
    int count1 = 0;
    for (int j = 0; j < header[3]/3; j++)
      {
      for (int k = 0; k < 3; k++)
        {
        faceId[k] = (int)data[(int)(header[0]*header[1]*header[2] + count1)];
        count1++;
        }
      //std::cout << "Loading (" << faceId[0] << "," << faceId[1] << "," << faceId[2] << ")" << std::endl;
      polys->InsertNextCell(3, faceId);
      }
    shape->SetPoints(points);
    points->Delete();
    shape->SetPolys(polys);
    polys->Delete();
    //std::cout << "Number of points in shape " << i << " is " << shape->GetNumberOfPoints() << std::endl;
    this->AddShape(shape);
    }
  this->m_WeightMode = false;
  this->m_RobustMode = false;
  this->m_StandardMode = true;

  delete [] faceId;
  delete [] data;
  delete [] header;
  itkDebugMacro(<< "Finished Reading SSM Data File");
  return true;
}

template <class TProfileSamplingPrecisionType>
bool
StatisticalShapeModel<TProfileSamplingPrecisionType>
::SaveModel(const char * filename)
{
  itkDebugMacro(<< "About to Save SSM with filename " << filename);
  // We want to save the current model to file
  // File format is header, shape (point) data, connectiveness (triangle) data
  int sz = this->m_Landmarks.size();
  //int sz = m_Landmarks->GetNumberOfItems();
  if (sz > 0)
    {
    std::ofstream fout(filename, std::ios::binary);
    if(!fout.is_open())
      {
      itkExceptionMacro(<< "Cannot create SSM3D File " << filename);
      return false;
      }

    int headerSize = 4;
    float * header = new float[headerSize];
    header[0] = sz;
    header[1] = 3;
    header[2] = this->GetShape(0)->GetPoints()->GetNumberOfPoints();
    header[3] = 3*this->GetShape(0)->GetPolys()->GetNumberOfCells();
    //std::cout << "Header is " << header[0] << " " << header[1] << " " << header[2] << std::endl;
    int sizeData = (int)(header[0]*header[1]*header[2] + header[3]);
    float *writer = new float [sizeData];

    int count = 0;
    double value[3];
    for (int i = 0; i < sz; i++)
      {
      vtkPolyData * shape = this->GetShape(i);
      for(int j = 0; j < shape->GetPoints()->GetNumberOfPoints(); j++)
        {
        shape->GetPoints()->GetPoint(j, value);
        writer[count] = value[0]; count++;
        writer[count] = value[1]; count++;
        writer[count] = value[2]; count++;
        }
      }
    int *fcarray = new int[3];
    // Added to remove 'false' valgrind errors about uninitialised data wrt below
    fcarray[0]= 0; fcarray[1]= 0; fcarray[2]= 0;
    //vtkIdType szv = 3;
    this->GetShape(0)->GetPolys()->InitTraversal();
    //std::cout << "Number of faces is " << m_Landmarks[0]->GetPolys()->GetNumberOfCells() << std::endl;
    vtkCellArray *cells = this->GetShape(0)->GetPolys();
    vtkTriangle *triangle;
    for(int i=0; i< cells->GetNumberOfCells(); i++)
      {
      triangle=(vtkTriangle *)this->GetShape(0)->GetCell(i);
      fcarray[0] = triangle->GetPointId(0);
      fcarray[1] = triangle->GetPointId(1);
      fcarray[2] = triangle->GetPointId(2);
      //cells->GetNextCell(szv,(vtkIdType *)fcarray);
      //std::cout << "Saving (" << fcarray[0] << "," << fcarray[1] << "," << fcarray[0] << ")"<< std::endl;
      writer[count] = fcarray[0]; count++;
      writer[count] = fcarray[1]; count++;
      writer[count] = fcarray[2]; count++;
      }
    delete [] fcarray;

    // Write data to file
    fout.write((char *)(header),headerSize*sizeof(float));
    fout.write((char *)(writer),sizeData*sizeof(float));
    fout.close();
    delete [] header;
    delete [] writer;
    itkDebugMacro(<< "Finished writing file " << filename);
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
StatisticalShapeModel<TProfileSamplingPrecisionType>
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
  if(this->m_PCA)
    this->m_PCA->Delete();
  this->m_PCA = vtkPCAAnalysisFilter::New(); //reset
}

}

#endif
