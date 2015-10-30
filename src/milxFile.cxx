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
#include "milxFile.h"
//ITK
#ifndef VTK_ONLY
  #include <itkAffineTransform.h>
  #include <itkVersorRigid3DTransform.h>
  #include <itkRigid3DTransform.h>
  #include <itkCenteredEuler3DTransform.h>
  #include <itkEuler3DTransform.h>
  #include <itkTransformFileWriter.h>
#endif
//VTK
#ifndef ITK_ONLY
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
  #include <vtkDelimitedTextReader.h>
  #include <vtkDelimitedTextWriter.h>
#endif

#include <zlib.h>

namespace milx
{

#ifndef VTK_ONLY
  bool File::CanReadImage(const std::string filename)
  {
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::ImageIOFactory::ReadMode);

    if(!imageIO)
      return false;
    else
      return true;
  }

  bool File::ReadImageInformation(const std::string filename, std::string &pixeltype, std::string &componentType, size_t &dimensions)
  {
    if(CanReadImage(filename))
    {
      itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::ImageIOFactory::ReadMode);
      imageIO->SetFileName(filename.c_str());
      try
      {
        imageIO->ReadImageInformation();
      }
      catch( itk::ExceptionObject & err )
      {
        std::cerr << "File Exception caught while reading header!" << std::endl;
        std::cerr << err << std::endl;
        return false;
      }

      pixeltype = imageIO->GetPixelTypeAsString(imageIO->GetPixelType());
      componentType = imageIO->GetComponentTypeAsString(imageIO->GetComponentType());
      dimensions = imageIO->GetNumberOfDimensions();
    }
    else
      return false;

    return true;
  }

  std::vector<std::string> File::GetSupportedImageFileExtensions()
  {
    std::vector<std::string> extensions;

    //Code taken from itk::ImageIOFactory->CreateImageIO() member in ITK 3.20
    std::list<itk::LightObject::Pointer> allobjects = itk::ObjectFactoryBase::CreateAllInstance("itkImageIOBase");
    for(std::list<itk::LightObject::Pointer>::iterator i = allobjects.begin(); i != allobjects.end(); ++i)
    {
      itk::ImageIOBase* imageIO = dynamic_cast<itk::ImageIOBase*>(i->GetPointer());
      if(imageIO)
      {
        for(size_t j = 0; j < imageIO->GetSupportedReadExtensions().size(); j ++)
          extensions.push_back(imageIO->GetSupportedReadExtensions()[j]);
      }
    }

    return extensions;
  }

  std::vector<std::string> File::GetDICOMSeriesUIDs(const std::string directoryPath, bool recursive)
  {
    std::vector<std::string> UIDs;
    typedef itk::GDCMSeriesFileNames GeneratorType;
    GeneratorType::Pointer nameGenerator = GeneratorType::New();
    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->AddSeriesRestriction("0008|0021");
    nameGenerator->SetDirectory(directoryPath.c_str());
    if(recursive)
        nameGenerator->RecursiveOn();
    nameGenerator->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  #if (ITK_VERSION_MAJOR > 3)
    try
    {
      nameGenerator->Update();
    }
    catch (itk::ExceptionObject &excp)
    {
      std::cerr << "File Exception caught while reading filenames in directory!" << std::endl;
      std::cerr << excp << std::endl;
      return UIDs;
    }
  #endif

    UIDs = nameGenerator->GetSeriesUIDs();

    return UIDs;
  }

  std::vector<std::string> File::GetDICOMSeriesFilenames(const std::string directoryPath, const std::string seriesName)
  {
    std::vector<std::string> filenames;
    typedef itk::GDCMSeriesFileNames GeneratorType;
    GeneratorType::Pointer nameGenerator = GeneratorType::New();
    nameGenerator->SetUseSeriesDetails(true);
    nameGenerator->AddSeriesRestriction("0008|0021");
    nameGenerator->SetDirectory(directoryPath.c_str());
    nameGenerator->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  #if (ITK_VERSION_MAJOR > 3)
    try
    {
      nameGenerator->Update();
    }
    catch (itk::ExceptionObject &excp)
    {
      std::cerr << "File Exception caught while reading filenames in directory!" << std::endl;
      std::cerr << excp << std::endl;
      return filenames;
    }
  #endif

    filenames = nameGenerator->GetFileNames(seriesName);

    return filenames;
  }
#endif

#ifndef ITK_ONLY
  #ifndef VTK_ONLY
  vtkMatrix4x4* File::OpenITKTransform(std::string filename)
  {
      vtkMatrix4x4 * matrix = vtkMatrix4x4::New();
      matrix->Identity();

      typedef itk::TransformFileReader        TransformReaderType;
      typedef TransformReaderType::TransformListType * TransformListType;
      typedef itk::AffineTransform< double, 3 > AffineTransformType;
      typedef AffineTransformType::Pointer AffineTransformPointer;
      typedef itk::VersorRigid3DTransform<double> VersorRigidTransformType;
      typedef VersorRigidTransformType::Pointer VersorRigidTransformPointer;
      typedef itk::Rigid3DTransform<double> RigidTransformType;
      typedef RigidTransformType::Pointer RigidTransformPointer;
      typedef itk::CenteredEuler3DTransform<double> CenteredEuler3DTransformType;
      typedef CenteredEuler3DTransformType::Pointer CenteredEuler3DTransformPointer;
      typedef itk::Euler3DTransform<double> Euler3DTransformType;
      typedef Euler3DTransformType::Pointer Euler3DTransformPointer;

      TransformReaderType::Pointer affineReader = TransformReaderType::New();
      affineReader->SetFileName(filename);
      affineReader->AddObserver(itk::ProgressEvent(), milx::ProgressUpdates);
      try
      {
          affineReader->Update();
      }
      catch(itk::ExceptionObject &ex)
      {
          PrintError("milxFile: Failed reading Transform " + std::string(ex.GetDescription()));
          return matrix;
      }
      TransformListType transforms = affineReader->GetTransformList();
      TransformReaderType::TransformListType::const_iterator tit = transforms->begin();
      AffineTransformPointer affineTransform;
      VersorRigidTransformPointer versorRigidTransform;
      RigidTransformPointer rigidTransform;
      CenteredEuler3DTransformPointer centeredEulerTransform;
      Euler3DTransformPointer euler3DTransform;

      if( !strcmp((*tit)->GetNameOfClass(),"AffineTransform") )
      {
          AffineTransformPointer affine_read = static_cast<AffineTransformType*>((*tit).GetPointer());
          affineTransform = dynamic_cast< AffineTransformType * >( affine_read.GetPointer() );
          if( affineTransform )
          {
              PrintInfo("Successfully Read Affine Transform file.");
          }
          else
          {
              PrintError("milxFile: Transform is not of Affine Transform Type. Returning Identity.");
              return matrix;
          }

          for (unsigned long i = 0; i < 3; i++)
          {
              matrix->SetElement(i, 0, affineTransform->GetMatrix()(i,0));
              matrix->SetElement(i, 1, affineTransform->GetMatrix()(i,1));
              matrix->SetElement(i, 2, affineTransform->GetMatrix()(i,2));
              matrix->SetElement(i, 3, affineTransform->GetOffset()[i]);
          }
      }
      else if( !strcmp((*tit)->GetNameOfClass(),"VersorRigid3DTransform") )
      {
          VersorRigidTransformPointer rigid_read = static_cast<VersorRigidTransformType*>((*tit).GetPointer());
          versorRigidTransform = dynamic_cast< VersorRigidTransformType * >( rigid_read.GetPointer() );
          if( versorRigidTransform )
          {
              PrintInfo("Successfully Read Versor Rigid Transform file.");
          }
          else
          {
              PrintError("milxFile: Transform Input is not of Versor Rigid Type. Returning Identity.");
              return matrix;
          }

          for (unsigned long i = 0; i < 3; i++)
          {
              matrix->SetElement(i, 0, versorRigidTransform->GetMatrix()(i,0));
              matrix->SetElement(i, 1, versorRigidTransform->GetMatrix()(i,1));
              matrix->SetElement(i, 2, versorRigidTransform->GetMatrix()(i,2));
              matrix->SetElement(i, 3, versorRigidTransform->GetOffset()[i]);
          }
      }
      else if( !strcmp((*tit)->GetNameOfClass(),"Rigid3DTransform") )
      {
        RigidTransformPointer rigid_read = static_cast<RigidTransformType*>((*tit).GetPointer());
        rigidTransform = dynamic_cast< RigidTransformType * >( rigid_read.GetPointer() );
        if( rigidTransform )
        {
            PrintInfo("Successfully Read Rigid Transform file.");
        }
        else
        {
            PrintError("milxFile: Transform Input is not of Rigid Type. Returning Identity.");
            return matrix;
        }

        for (unsigned long i = 0; i < 3; i++)
        {
            matrix->SetElement(i, 0, rigidTransform->GetMatrix()(i,0));
            matrix->SetElement(i, 1, rigidTransform->GetMatrix()(i,1));
            matrix->SetElement(i, 2, rigidTransform->GetMatrix()(i,2));
            matrix->SetElement(i, 3, rigidTransform->GetOffset()[i]);
        }
      }
      else if( !strcmp((*tit)->GetNameOfClass(),"CenteredEuler3DTransform") )
      {
        CenteredEuler3DTransformPointer rigid_read = static_cast<CenteredEuler3DTransformType*>((*tit).GetPointer());
        centeredEulerTransform = dynamic_cast< CenteredEuler3DTransformType * >( rigid_read.GetPointer() );
        if( centeredEulerTransform )
        {
            PrintInfo("Successfully Read Centered Euler 3D Transform file.");
        }
        else
        {
            PrintError("milxFile: Transform Input is not of Centered Euler 3D Type. Returning Identity.");
            return matrix;
        }

        for (unsigned long i = 0; i < 3; i++)
        {
            matrix->SetElement(i, 0, centeredEulerTransform->GetMatrix()(i,0));
            matrix->SetElement(i, 1, centeredEulerTransform->GetMatrix()(i,1));
            matrix->SetElement(i, 2, centeredEulerTransform->GetMatrix()(i,2));
            matrix->SetElement(i, 3, centeredEulerTransform->GetOffset()[i]);
        }
      }
      else if( !strcmp((*tit)->GetNameOfClass(),"Euler3DTransform") )
      {
        Euler3DTransformPointer rigid_read = static_cast<Euler3DTransformType*>((*tit).GetPointer());
        euler3DTransform = dynamic_cast< Euler3DTransformType * >( rigid_read.GetPointer() );
        if( euler3DTransform )
        {
            PrintInfo("Successfully Read Euler 3D Transform file.");
        }
        else
        {
            PrintError("milxFile: Transform Input is not of Euler 3D Type. Returning Identity.");
            return matrix;
        }

        for (unsigned long i = 0; i < 3; i++)
        {
            matrix->SetElement(i, 0, euler3DTransform->GetMatrix()(i,0));
            matrix->SetElement(i, 1, euler3DTransform->GetMatrix()(i,1));
            matrix->SetElement(i, 2, euler3DTransform->GetMatrix()(i,2));
            matrix->SetElement(i, 3, euler3DTransform->GetOffset()[i]);
        }
      }
      else
      {
          PrintError("milxFile: Transform Input is not of known type");
      }

      return matrix;
  }

  vtkMatrix4x4* File::OpenVTKTransform(std::string filename)
  {
      vtkMatrix4x4 * matrix = vtkMatrix4x4::New();

      char * buffer = new char[512];
      gzFile fin = ::gzopen( filename.c_str(), "rb" );
      if(fin == NULL)
      {
          cerr << "Cannot read " << filename << endl;
          delete [] buffer;
      }
      buffer = ::gzgets (fin, buffer, 512);
      if(buffer[0] != '(')
      {
          cerr << "File is not a transform file" << endl;
      }
      buffer = ::gzgets (fin, buffer, 512);
      if((buffer[0] != 0) && (buffer[1] != '8'))
      {
          // Failed Magic number test ;)
          cerr << buffer << endl;
          cerr << "File is not a transform file" << endl;
      }
      char str1 [128];
      char str2 [128];
      char str3 [128];
      char str4 [128];
      for (unsigned long i = 0; i < 4; i++)
      {
          //fin >> value;
          buffer = ::gzgets (fin, buffer, 512);

          //sscanf(buffer,"%f %f %f", &x, &y, &z);
          sscanf(buffer,"%s %s %s %s", &str1[0], &str2[0], &str3[0], &str4[0]);
          matrix->SetElement(i, 0, atof(str1));
          matrix->SetElement(i, 1, atof(str2));
          matrix->SetElement(i, 2, atof(str3));
          matrix->SetElement(i, 3, atof(str4));
          std::cout << atof(str1) << " " << atof(str2) << " " << atof(str3) << " " << atof(str4) << std::endl;
      }
      buffer = ::gzgets (fin, buffer, 512);
      if(buffer[0] != ')')
      {
          cerr << buffer << endl;
          cerr << "File is not a transform file" << endl;
      }
      delete [] buffer;
      gzclose(fin);
      return matrix;
  }

  void File::SaveITKTransform(std::string filename, const vtkMatrix4x4 * matrix)
  {
    // Step 1.: Convert VTK matrix to an ITK transform
    typedef itk::AffineTransform<double,3> TransformType;
    TransformType::Pointer transform = TransformType::New();
    TransformType::MatrixType itkmat;

    for (size_t i=0; i<3; i++) {
      for (size_t j=0; j<3; j++) {
        itkmat(i,j) = matrix->GetElement(i,j);
      }
    }
    transform->SetMatrix(itkmat);

    TransformType::OutputVectorType trans;
    trans[0] = matrix->GetElement(0,3);
    trans[1] = matrix->GetElement(1,3);
    trans[2] = matrix->GetElement(2,3);
    transform->SetTranslation(trans);

    // Step 2.: Save
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter->SetFileName(filename);
    transformWriter->SetInput(transform);
    try {
      transformWriter->Update();
    } catch( itk::ExceptionObject & err) {
      std::cerr << "Exception occurred while writing transform file:\n    " << filename << std::endl;
      std::cerr << err << std::endl;
      std::cout << err.what() << std::endl;
    }
  }

  /*The file syntax is like, eg:
      (
      O8
      1.000000 0.000000  0.000000  22.500000
      0.000000 1.000000  0.000000  22.500000
      0.000000 0.000000  1.000000  22.500000
      0.000000 0.000000  0.000000   1.000000
      )
      Note: 08 is a magic number.
  */
  void File::SaveVTKTransform(std::string filename, const vtkMatrix4x4 * matrix)
  {
    std::string ext = milx::File::GetFileExtension(filename);

    if (ext == "gz")
    {
      gzFile fout = ::gzopen( filename.c_str(), "wb" );

      if(fout == NULL)
      {
        cerr << "Cannot open for writing: " << filename << endl;
      }

      ::gzputs(fout, "(\n08\n");

      for (unsigned long i = 0; i < 4; i++)
      {
        ::gzprintf(fout, "%- 24.18g %- 24.18g %- 24.18g %- 24.18g\n", matrix->GetElement(i, 0), matrix->GetElement(i, 1), matrix->GetElement(i, 2), matrix->GetElement(i, 3));
      }

      ::gzputs(fout, ")\n");
      ::gzclose(fout);
    }
    else
    {
      FILE * fout = ::fopen( filename.c_str(), "wb" );

      if(fout == NULL)
      {
        cerr << "Cannot open for writing: " << filename << endl;
      }

      ::fputs("(\n08\n", fout);

      for (unsigned long i = 0; i < 4; i++)
      {
        ::fprintf(fout, "%- 24.18g %- 24.18g %- 24.18g %- 24.18g\n", matrix->GetElement(i, 0), matrix->GetElement(i, 1), matrix->GetElement(i, 2), matrix->GetElement(i, 3));
      }

      ::fputs(")\n", fout);
      ::fclose(fout);
    }


  }
  #endif // VTK_ONLY

bool File::OpenModel(const std::string filename, vtkSmartPointer<vtkPolyData> &data)
{
  bool legacy = false, wavefront = false, stanford = false, stereoLith = false;
  std::string extension = GetFileExtension(filename);

  if (extension == "vtk")
    legacy = true; //Load legacy VTK file
  else if (extension == "obj")
    wavefront = true;
  else if (extension == "ply")
    stanford = true;
  else if(extension == "stl")
    stereoLith = true;

  if(data != NULL)
    PrintError("milxFile: PolyData pointer is not NULL. May get a memory leak.");

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

  PrintError("Reader Encountered the following error.");
  PrintError(errorObserver->GetMessage());
  return false;
}

bool File::OpenModelCollection(std::vector<std::string> filenames, vtkSmartPointer<vtkPolyDataCollection> &collection)
{
  collection = vtkSmartPointer<vtkPolyDataCollection>::New();

  for(std::vector<std::string>::iterator name = filenames.begin(); name != filenames.end(); name ++)
  {
    vtkSmartPointer<vtkPolyData> surface;

    if( !File::OpenModel(*name, surface) ) //Error printed inside
      continue;

    collection->AddItem( surface );
  }

  return true;
}

bool File::SaveModel(const std::string filename, vtkSmartPointer<vtkPolyData> data, const bool binary)
{
  bool legacy = false, stanford = false, stereoLith = false;
  std::string extension = GetFileExtension(filename);

  if (extension == "vtk")
    legacy = true; //Save legacy VTK file
  else if (extension == "ply")
    stanford = true; //Save Stanford Poly file
  else if(extension == "stl")
    stereoLith = true; //STL files

  if(data == NULL)
  {
    PrintError("PolyData pointer is NULL. Not saving.");
    return false;
  }

  if (legacy)
  {
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
      writer->SetFileName(filename.c_str());
    #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(data);
    #else
      writer->SetInputData(data);
    #endif
      writer->Write();
  }
  else if (stanford)
  {
    vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
      writer->SetFileName(filename.c_str());
    #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(data);
    #else
      writer->SetInputData(data);
    #endif
      writer->Write();
  }
  else if(stereoLith)
  {
    vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
      writer->SetFileName(filename.c_str());
    #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(data);
    #else
      writer->SetInputData(data);
    #endif
      writer->Write();
  }
  else
  {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      writer->SetFileName(filename.c_str());
    #if VTK_MAJOR_VERSION <= 5
      writer->SetInput(data);
    #else
      writer->SetInputData(data);
    #endif
      if (binary)
        writer->SetDataModeToBinary();
      writer->SetCompressorTypeToZLib();
      writer->Write();
  }

  return true;
}

bool File::SaveModelCollection(std::vector<std::string> filenames, vtkSmartPointer<vtkPolyDataCollection> collection, const bool binary)
{
  collection->InitTraversal();
  for(std::vector<std::string>::iterator name = filenames.begin(); name != filenames.end(); name ++)
  {
    if( !File::SaveModel(*name, collection->GetNextItem(), binary) ) //Error printed inside
      continue;
  }

  return true;
}

bool File::SaveTransform(const std::string filename, vtkSmartPointer<vtkTransform> data, const bool ITK)
{
  if(data == NULL)
  {
    PrintError("Transform pointer is NULL. Not saving.");
    return false;
  }

  if (ITK) {
    SaveITKTransform(filename, data->GetMatrix());
  } else {
    SaveVTKTransform(filename, data->GetMatrix());
  }

  return true;
}

bool File::SaveTransformCollection(std::vector<std::string> filenames, vtkSmartPointer<vtkTransformCollection> collection, const bool ITK)
{
  collection->InitTraversal();
  for(std::vector<std::string>::iterator name = filenames.begin(); name != filenames.end(); name ++)
  {
    if( !File::SaveTransform(*name, collection->GetNextItem(), ITK) ) //Error printed inside
      continue;
  }

  return true;
}

bool File::OpenDelimitedText(const std::string filename, vtkSmartPointer<vtkTable> &data, const std::string delimiter)
{
  PrintInfo("Reading file: " + filename);
  vtkSmartPointer<vtkDelimitedTextReader> csv_vert_source = vtkSmartPointer<vtkDelimitedTextReader>::New();
    csv_vert_source->SetFieldDelimiterCharacters(delimiter.c_str());
    csv_vert_source->DetectNumericColumnsOn();
    csv_vert_source->SetHaveHeaders(true);
    csv_vert_source->SetFileName(filename.c_str());
    csv_vert_source->Update();
    PrintInfo("Has " + NumberToString(csv_vert_source->GetOutput()->GetNumberOfColumns()) + " columns");

  data = csv_vert_source->GetOutput();

  return true;
}

bool File::SaveDelimitedText(const std::string filename, const vtkSmartPointer<vtkTable> data, const std::string delimiter)
{
  PrintInfo("Writing file: " + filename);
  vtkSmartPointer<vtkDelimitedTextWriter> csv_vert_source = vtkSmartPointer<vtkDelimitedTextWriter>::New();
    csv_vert_source->SetFieldDelimiter(delimiter.c_str());
    csv_vert_source->SetFileName(filename.c_str());
  #if VTK_MAJOR_VERSION <= 5
    csv_vert_source->SetInput(data);
  #else
    csv_vert_source->SetInputData(data);
  #endif
    csv_vert_source->Write();

  return true;
}

bool File::SaveCamera(const std::string filename, const vtkSmartPointer<vtkCamera> camera)
{
  if(!camera || filename.empty())
    return false;

  ofstream file(filename.c_str(), ios::out);

  if(file.fail())
    return false;

  file << camera->GetClippingRange()[0] << " "
    << camera->GetClippingRange()[1] << " "

    << camera->GetFocalPoint()[0] << " "
    << camera->GetFocalPoint()[1] << " "
    << camera->GetFocalPoint()[2] << " "

    << camera->GetPosition()[0] << " "
    << camera->GetPosition()[1] << " "
    << camera->GetPosition()[2] << " "

    << camera->GetViewAngle() << " "

    << camera->GetViewUp()[0] << " "
    << camera->GetViewUp()[1] << " "
    << camera->GetViewUp()[2] << " "
    << std::endl;
  file.close();

  return true;
}

bool File::LoadCamera(const std::string filename, vtkSmartPointer<vtkCamera> &camera)
{
  if(!camera || filename.empty())
    return false;

  double clippingRange[2];
	double focal[3];
	double pos[3];
	double viewAngle;
	double up[3];

	ifstream file (filename.c_str());

	if(file.fail())
    return false;

	file
    >> clippingRange[0]
    >> clippingRange[1]

    >> focal[0]
    >> focal[1]
    >> focal[2]

    >> pos[0]
    >> pos[1]
    >> pos[2]

    >> viewAngle

    >> up[0]
    >> up[1]
    >> up[2];
	file.close();

	// set the camera details
	camera->SetClippingRange(clippingRange);
	camera->SetFocalPoint(focal);
	camera->SetPosition(pos);
	camera->SetViewAngle(viewAngle);
	camera->SetViewUp(up);

	return true;
}
#endif //ITK_ONLY

} //end milx namespace
