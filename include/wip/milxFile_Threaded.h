/*=========================================================================
  Copyright: (c) CSIRO, Australia

  This software is protected by international copyright laws.
  Any unauthorised copying, distribution or reverse engineering is prohibited.

  Licence:
  All rights in this Software are reserved to CSIRO. You are only permitted
  to have this Software in your possession and to make use of it if you have
  agreed to a Software License with CSIRO.
=========================================================================*/

#ifndef __MILXFILE_H
#define __MILXFILE_H

//STL
#include <string>
#include <algorithm>
#include <vector>
//ITK
#include <itkImage.h>
#include <itkMultiThreader.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
//VTK
#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkImageData.h>
#include <vtkImageFlip.h>
#include <vtkPolyData.h>
#include <vtkPolyDataCollection.h>
//MILX
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
//SMILI
#include <milxGlobal.h>

typedef std::string string; //Don't open std namespace

/**
  \class vtkErrorWarning
  \brief Object for intercepting errors thrown by VTK based readers.
*/
class VTK_COMMON_EXPORT vtkErrorWarning : public vtkCommand
{
public:
  /**
      \fn vtkErrorWarning::New()
      \brief Standard VTK Object Factory
  */
  static vtkErrorWarning* New()
  {
    return new vtkErrorWarning();
  }

  vtkTypeMacro(vtkErrorWarning, vtkCommand);

  inline const char* GetMessage()
  {
    return m_Message.c_str();
  }

  inline bool HasFailed()
  {
    return m_ErrorEncountered;
  }

  inline bool ReportsFailure()
  {
    return m_ErrorEncountered;
  }

  void Initialize()
  {
    m_ErrorEncountered = m_WarningEncountered = false;
  }

protected:
  vtkErrorWarning();
  ~vtkErrorWarning();

  bool m_ErrorEncountered;
  bool m_WarningEncountered;
  string m_Message;

private:
  void Execute(vtkObject *caller, unsigned long observedType, void* message);
};

namespace milx
{

/**
  \class File
  \brief A general file IO class for reading writing images and models.
  
  Loading a model:
  \code
  vtkSmartPointer<vtkPolyData> surface;
  std::string name = getFilenameFromSomewhere();

  if( !File::OpenModel(name, surface) ) //Error printed inside
    return false;
  \endcode
  Loading models:
  \code
  std::vector<std::string> filenames = getFilenamesFromSomewhere();
  vtkSmartPointer<vtkPolyDataCollection> collection;
  if( !milx::File::OpenModelCollection(filenames, collection) ) //Error printed inside
    exit(EXIT_FAILURE);
  \endcode
*/
//template<class Type = float>
class SMILI_EXPORT File
{
public:
  /*!
    \fn File::File()
    \brief Standard constructor
	*/
  File() {};
  /*!
    \fn File::~File()
    \brief Standard Destructor
	*/
  ~File() {};

public:
  //Images
  /*!
    \fn File::OpenImage(const string filename, typename itk::Image<Type, Dimension>::Pointer &data)
    \brief Opens an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI etc.

    Returns true if successful. data is allocated within this member, so pass a NULL pointer.
    Image is also NOT flipped, consider using overloaded OpenImage() with VTK image data which is flipped.
  */
  template<class Type, unsigned int Dimension>
  bool OpenImage(const string filename, typename itk::Image<Type, Dimension>::Pointer &data);
  /*!
    \fn File::OpenImage(const string filename, vtkSmartPointer<vtkImageData> &data)
    \brief Opens an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI, INR etc. Overloaded for VTK input/output.

    ITK Reader is used and then transparently converted to a VTK image.
    Returns true if successful. Image is also flipped since the ITK reader orientation is different to VTK images.
  */
  template<class Type, unsigned int Dimension>
  bool OpenImage(const string filename, vtkSmartPointer<vtkImageData> &data);
  /*!
    \fn File::SaveImage(const string filename, typename itk::Image<Type, Dimension>::Pointer data)
    \brief Saves an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI etc.

    Returns true if successful.
    Image is also NOT flipped, consider using overloaded OpenImage() with VTK image data which is flipped.
  */
  template<class Type, unsigned int Dimension>
  bool SaveImage(const string filename, typename itk::Image<Type, Dimension>::Pointer data);
  /*!
    \fn File::SaveImage(const string filename, vtkSmartPointer<vtkImageData> data)
    \brief Saves an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI, INR etc. Overloaded for VTK input/output.

    ITK Writer is used after transparently converting from a VTK image.
    Returns true if successful. Image is also flipped since the ITK reader orientation is different to VTK images.
  */
  template<class Type, unsigned int Dimension>
  bool SaveImage(const string filename, vtkSmartPointer<vtkImageData> data);

  //Surfaces/Models
  /*!
    \fn File::OpenModel(const string filename, vtkSmartPointer<vtkPolyData> &data)
    \brief Opens a model file, which can be a VTK XML, Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk), a Polygonal File (*.ply) or a Object file (*.obj).

    Returns true if successful. PolyData is allocated within this member, so pass a NULL pointer. Hence, the pass is by reference.
  */
  static bool OpenModel(const string filename, vtkSmartPointer<vtkPolyData> &data);
  /*!
    \fn File::OpenModelCollection(std::vector<string> filenames, vtkSmartPointer<vtkPolyDataCollection> &collection)
    \brief Opens model files, which can be a VTK XML, Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk), a Polygonal File (*.ply) or a Object file (*.obj) into a collection.

    Returns true if successful. PolyData and PolyDataCollection are allocated within this member, so pass a NULL pointer. Hence, the pass is by reference.
  */
  static bool OpenModelCollection(std::vector<string> filenames, vtkSmartPointer<vtkPolyDataCollection> &collection);
  /*!
    \fn DGVFile::SaveModel(const string filename, vtkSmartPointer<vtkPolyData> data, const bool binary = false)
    \brief Saves a model as a file, which can be a VTK XML or Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk) or a Polygonal File (*.ply).

    Returns true if successful.
    */
  static bool SaveModel(const string filename, vtkSmartPointer<vtkPolyData> data, const bool binary = false);
  /*!
    \fn File::SaveModelCollection(std::vector<string> filenames, vtkSmartPointer<vtkPolyDataCollection> collection)
    \brief Saves model files, which can be a VTK XML, Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk), a Polygonal File (*.ply) or a Object file (*.obj) from a collection.

    Returns true if successful. 
  */
  static bool SaveModelCollection(std::vector<string> filenames, vtkSmartPointer<vtkPolyDataCollection> collection);

  /*!
    \fn File::StringToLowerCase(string filename)
    \brief Returns a lower case version of the string.
  */
  static inline string StringToLowerCase(string filename)
  {
    std::transform(filename.begin(), filename.end(), filename.begin(), ::tolower);
    return filename;
  }
  /*!
    \fn File::StripFileExtension(const string filename)
    \brief Returns filename with last extension removed.
  
    "/foo/bar/baz.txt" --> "txt"
  */
  static inline string StripFileExtension(const string &filename)
  {
    return filename.substr(0, filename.find_last_of("."));
  }
  /*!
    \fn File::ExtractFilename(const string &filename, char delimiter = '/')
    \brief Returns filename (with extension) with paths removed.
  
    e.g. "/foo/bar/baz.txt" --> "baz.txt"
  */
  static inline string ExtractFilename(const string &filename, char delimiter = '/')
  {
    return filename.substr( filename.find_last_of( delimiter ) + 1 );
  }
  /*!
    \fn File::ExtractPath(const string &filename, char delimiter = '/')
    \brief Returns the path of the filename with the basename removed.
  
    "/foo/bar/baz.txt" --> "/foo/bar/"
  */
  static inline string ExtractPath(const string &filename, char delimiter = '/')
  {
    return filename.substr( 0, filename.find_last_of( delimiter ) + 1 );
  }
  /*!
    \fn File::GetBaseName(string filename)
    \brief Returns filename with paths and extension removed.
  
    e.g. "/foo/bar/baz.txt" --> "baz"
  */
  static inline string GetBaseName(const string &filename)
  {
    return StripFileExtension( ExtractFilename(filename) );
  }
  /*!
    \fn File::GetFileExtension(const string &filename)
    \brief Returns the file extension (in lower case) of the given filename

    This uses the last occurance of the period to get the right most substring.
  */
  static inline string GetFileExtension(const string &filename)
  {
    if (filename.find_last_of(".") != string::npos)
      return StringToLowerCase( filename.substr(filename.find_last_of(".")+1) );
    return "";
  }

protected:
  /**
  * \struct ThreadStruct
  * \brief Thread Struct for running IO classes in threads
  */
  struct ThreadStruct
  {
    typename itk::ProcessObject::Pointer Filter;
  };
  
  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to ThreadedGenerateData(). */
  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );
  
private:

};

// Callback routine used by the threading library. This routine just calls
// the ThreadedGenerateData method after setting the correct region for this
// thread. 
ITK_THREAD_RETURN_TYPE File::ThreaderCallback( void *arg )
{
  ThreadStruct *str;

  str = (ThreadStruct *)(((itk::MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  try
  {
    str->Filter->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    PrintError("File Exception caught!");
    err.Print(cerr);
  }
  
  return ITK_THREAD_RETURN_VALUE;
}

template<class Type, unsigned int Dimension>
bool File::OpenImage(const string filename, typename itk::Image<Type, Dimension>::Pointer &data)
{
  typedef itk::ImageFileReader< itk::Image<Type, Dimension> > ImageReader;
  
  /** Internal structure used for passing image data into the threading library */
  ThreadStruct threadData;
  
  typename ImageReader::Pointer reader = ImageReader::New();
    reader->SetFileName(filename.c_str());
    threadData.Filter = reader;
    
  //Threading
  typedef itk::MultiThreader ThreaderType;
  
  ThreaderType::Pointer threader = ThreaderType::New();
    threader->SetNumberOfThreads( 1 );
    threader->SetSingleMethod( ThreaderCallback, &threadData );
    threader->SingleMethodExecute();  
    
  data = reader->GetOutput();
  
  return true;
}

template<class Type, unsigned int Dimension>
bool File::OpenImage(const string filename, vtkSmartPointer<vtkImageData> &data)
{
  typename itk::Image<Type, Dimension>::Pointer tmpImage;
  
  if( !OpenImage(filename, tmpImage) )
    return false;
  
  typedef itk::ImageToVTKImageFilter< itk::Image<Type, Dimension> > ConnectorType;
  ///Export to VTK
  typename ConnectorType::Pointer connector = ConnectorType::New(); //!< ITK to VTK image
    connector->SetInput(tmpImage);
    connector->Update();
  
  ///Flip - ITK and VTK do not have the same orientation
  vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
    imageReorient->SetInput(connector->GetOutput());
    imageReorient->SetFilteredAxis(1);
    imageReorient->FlipAboutOriginOn();
    imageReorient->Update();
  
  data = imageReorient->GetOutput();
  
  return true;
}

template<class Type, unsigned int Dimension>
bool File::SaveImage(const string filename, typename itk::Image<Type, Dimension>::Pointer data)
{
  typedef itk::ImageFileWriter< itk::Image<Type, Dimension> > ImageWriter;
  
  /** Internal structure used for passing image data into the threading library */
  ThreadStruct threadData;
  
  typename ImageWriter::Pointer writer = ImageWriter::New();
    writer->SetFileName(filename.c_str());
    threadData.Filter = writer;
    
  //Threading
  typedef itk::MultiThreader ThreaderType;
  
  ThreaderType::Pointer threader = ThreaderType::New();
    threader->SetNumberOfThreads( 1 );
    threader->SetSingleMethod( ThreaderCallback, &threadData );
    threader->SingleMethodExecute();  
    
  return true;
}

template<class Type, unsigned int Dimension>
bool File::SaveImage(const string filename, vtkSmartPointer<vtkImageData> data)
{
  ///Flip - ITK and VTK do not have the same orientation
  vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
    imageReorient->SetInput(data);
    imageReorient->SetFilteredAxis(1);
    imageReorient->FlipAboutOriginOn();
    imageReorient->Update();
  
  typedef itk::VTKImageToImageFilter< itk::Image<Type, Dimension> > ConnectorType;
  ///Export to VTK
  typename ConnectorType::Pointer connector = ConnectorType::New(); //!< ITK to VTK image
    connector->SetInput(imageReorient->GetOutput());
    connector->Update();
  
  if( !SaveImage(filename, connector->GetOutput()) )
    return false;
  
  return true;
}

} //end namespace milx

#endif //__MILXFILE_H
