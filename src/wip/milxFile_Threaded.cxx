/*=========================================================================
  Program: milxSMILI
  Module: milxModel
  Author: Shekhar Chandra
  Language: C++
  Created: 18 Apr 2011 12:03:00 EST

  Copyright: (c) CSIRO, Australia

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.

  BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/
#include "milxFile.h"
//VTK
#include <vtkOBJReader.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>

namespace milx
{
  
bool File::OpenModel(const string filename, vtkSmartPointer<vtkPolyData> &data)
{
  bool legacy = false, wavefront = false, stanford = false, stereoLith = false;
  string extension = GetFileExtension(filename);

  if (extension == "vtk")
    legacy = true; //Load legacy VTK file
  else if (extension == "obj")
    wavefront = true;
  else if (extension == "ply")
    stanford = true;
  else if(extension == "stl")
    stereoLith = true;
  
  if(data != NULL)
    std::cerr << "PolyData pointer is not NULL. May get a memory leak." << endl;

  vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();
  if (legacy)
  {
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    reader->Update();

    if (!errorObserver->ReportsFailure())
    {
      data = reader->GetOutput();
      return true;
    }
  }
  else if (wavefront)
  {
    vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(filename.c_str());
    reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    reader->Update();

    if (!errorObserver->ReportsFailure())
    {
      data = reader->GetOutput();
      return true;
    }
  }
  else if (stanford)
  {
    vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(filename.c_str());
    reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    reader->Update();

    if (!errorObserver->ReportsFailure())
    {
      data = reader->GetOutput();
      return true;
    }
  }
  else if(stereoLith)
  {
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(filename.c_str());
    reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    reader->Update();

    if (!errorObserver->ReportsFailure())
    {
      data = reader->GetOutput();
      return true;
    }
  }
  else
  {
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    reader->Update();
    
    if (!errorObserver->ReportsFailure())
    {
      data = reader->GetOutput();
      return true;
    }
  }
  
  cerr << "Reader Encountered the following error." << endl;
  cerr << errorObserver->GetMessage() << endl;
  return false;
}

bool File::OpenModelCollection(std::vector<string> filenames, vtkSmartPointer<vtkPolyDataCollection> &collection)
{
  collection = vtkSmartPointer<vtkPolyDataCollection>::New();
  
  for(std::vector<string>::iterator name = filenames.begin(); name != filenames.end(); name ++)
  {
    vtkSmartPointer<vtkPolyData> surface;

    if( !File::OpenModel(*name, surface) ) //Error printed inside
      return false;
    
    collection->AddItem( surface );
  }
  
  return true;
}

bool File::SaveModel(const string filename, vtkSmartPointer<vtkPolyData> data, const bool binary)
{
  bool legacy = false, stanford = false, stereoLith = false;
  string extension = GetFileExtension(filename);

  if (extension == "vtk")
    legacy = true; //Save legacy VTK file
  else if (extension == "ply")
    stanford = true; //Save Stanford Poly file
  else if(extension == "stl")
    stereoLith = true; //STL files
  
  if(data == NULL)
  {
    std::cerr << "PolyData pointer is NULL. Not saving." << endl;
    return false;
  }

  if (legacy)
  {
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(data);
    writer->Write();
  }
  else if (stanford)
  {
    vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(data);
    writer->Write();
  }
  else if(stereoLith)
  {
    vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(data);
    writer->Write();
  }
  else
  {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(data);
    if (binary)
      writer->SetDataModeToBinary();
    writer->Write();
  }

  return true;
}

bool File::SaveModelCollection(std::vector<string> filenames, vtkSmartPointer<vtkPolyDataCollection> collection)
{
  collection->InitTraversal();
  for(std::vector<string>::iterator name = filenames.begin(); name != filenames.end(); name ++)
  {
    if( !File::SaveModel(*name, collection->GetNextItem()) ) //Error printed inside
      return false;
  }
  
  return true;
}

} //end milx namespace

vtkErrorWarning::vtkErrorWarning()
{
  Initialize();
}

vtkErrorWarning::~vtkErrorWarning()
{
  
}

void vtkErrorWarning::Execute(vtkObject *caller, unsigned long observedType, void* message)
{
    m_Message = (const char*) message;

    if(observedType == vtkCommand::ErrorEvent)
        m_ErrorEncountered = true;
    if(observedType == vtkCommand::WarningEvent)
        m_WarningEncountered = true;
}
