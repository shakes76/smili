/*=========================================================================
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
#ifndef __MILXFILE_H
#define __MILXFILE_H

//STL
#include <utility>
#include <string>
#include <algorithm>
#include <vector>
#include <sys/stat.h> //exists check
//ITK
#ifndef VTK_ONLY
  #include <itkImage.h>
  #include <itkImageFileReader.h>
  #include <itkImageSeriesReader.h>
  #include <itkImageFileWriter.h>
  #include <itkTransformFileReader.h>
  #include <itkTransform.h>
  #include <itkTransformFactoryBase.h>
  #include <itkGDCMImageIO.h>
  #include <itkGDCMSeriesFileNames.h>
  #include <itkOrientImageFilter.h> //for dicom orientation
  #include <itkExtractImageFilter.h>
  #if (ITK_VERSION_MAJOR > 3)
    #include <itkComposeImageFilter.h>
  #endif
#endif
//VTK
#ifndef ITK_ONLY
  #include <vtkSmartPointer.h>
  #include <vtkCommand.h>
  #include <vtkMatrix4x4.h>
  #include <vtkImageData.h>
  #include <vtkImageFlip.h>
  #include <vtkPolyData.h>
  #include <vtkPolyDataCollection.h>
  #include <vtkTransformCollection.h>
  #include <vtkTable.h>
  #include <vtkCamera.h>
  #include <vtkXMLImageDataReader.h>
  #include <vtkXMLImageDataWriter.h>
  //MILX
  #ifndef VTK_ONLY
    #include "itkImageToVTKImageFilter.h"
    #include "itkVTKImageToImageFilter.h"
  #endif
#endif
//SMILI
#include <milxGlobal.h>

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
  virtual ~File() {};

public:
  //Images
#ifndef VTK_ONLY
  /*!
    \fn File::CanReadImage(const std::string filename)
    \brief Checks if the image can be read by supported libraries.

    Currently does ITK only.
  */
  static bool CanReadImage(const std::string filename);
  /*!
    \fn File::ReadImageInformation(const std::string filename, std::string &pixeltype, std::string &componentType, size_t &dimensions)
    \brief Reads just the header of an image file (without reading the image data), and writes the info into the arguments passed by reference such as type etc.

    \code
    //Check type etc. of medical image
    std::string pixelType, componentType;
    if(!milx::File::ReadImageInformation(filename, pixelType, componentType, dataDimensions))
    {
        cerr << "Failed reading header of image. File may not be an image. Exiting" << endl;
        return false;
    }
    \endcode
  */
  static bool ReadImageInformation(const std::string filename, std::string &pixeltype, std::string &componentType, size_t &dimensions);
  /*!
    \fn File::GetSupportedImageFileExtensions()
    \brief Determines the supported image file format as file extensions from internal libraries and returns a list of strings for each format.

    Currently does ITK only.
  */
  static std::vector<std::string> GetSupportedImageFileExtensions();
  /*!
    \fn File::OpenImage(const std::string filename, typename itk::SmartPointer<TImage> &data)
    \brief Opens an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI etc.

    Returns true if successful. data is allocated within this member, so pass a NULL pointer.
    Image is also NOT flipped, consider using overloaded OpenImage() with VTK image data which is flipped.
  */
  template<class TImage>
  static bool OpenImage(const std::string filename, typename itk::SmartPointer<TImage> &data);
#if (ITK_VERSION_MAJOR > 3)
  /*!
    \fn File::OpenAsVectorImage(const std::string filename, typename itk::SmartPointer< itk::VectorImage< typename TImage::InternalPixelType, typename TImage::ImageDimension-1> > &data)
    \brief Force open an N-D image file as a (N-1)-D vector image, i.e. the higher dimension is loaded into a variable length vector at each pixel.

    The image template parameter should be of N-D image and the dimension parameter should be N-1.
    \code
    typedef itk::VectorImage<float, 3> vectorImageType;
    vectorImageType::Pointer vecImg;

    typedef itk::Image<floatPixelType, 4> float4DImageType;
    if( !milx::File::OpenAsVectorImage<float4DImageType, 3>(filename.toStdString(), vecImg) )
        return false;
    \endcode

    This member requires ITK 4 or above.
    Returns true if successful. data is allocated within this member, so pass a NULL pointer.
    Image is also NOT flipped, consider using overloaded OpenImage() with VTK image data which is flipped.
  */
  template<class TImage, size_t TDim>
  static bool OpenAsVectorImage(const std::string filename, typename itk::SmartPointer< itk::VectorImage< typename TImage::InternalPixelType, TDim> > &data);
#endif
  /*!
    \fn File::OpenImages(std::vector<string> filenames, std::vector< typename itk::SmartPointer<TImage> > images)
    \brief Opens a number of images, which are any of the following: JPEG, PNG, DICOM, TIFF, NIFTI etc.

    Returns true if successful. images vector is filled with images. Image is also NOT flipped for VTK.
    The filenames vector is reassigned to the names that were actually read in (indices corresponding to the vector of images)
  */
  template<class TImage>
  static bool OpenImages(std::vector<std::string> &filenames, std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
    \fn File::SaveImage(const std::string filename, typename itk::SmartPointer<TImage> data, itk::ImageIOBase *io = NULL)
    \brief Saves an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI etc.

    Returns true if successful.
    Image is also NOT flipped, consider using overloaded OpenImage() with VTK image data which is flipped.
  */
  template<class TImage>
  static bool SaveImage(const std::string filename, typename itk::SmartPointer<TImage> data, itk::ImageIOBase *io = NULL);
  /*!
    \fn File::SaveImages(std::vector<std::string> &filenames, const std::vector< typename itk::SmartPointer<TImage> > images)
    \brief Saves a number of images, which are any of the following: JPEG, PNG, DICOM, TIFF, NIFTI etc.

    Returns true if successful. images vector is assumed to be filled with images. Image is also NOT flipped for VTK.
    The filenames vector is reassigned to the names that were actually written (indices corresponding to the vector of images)
  */
  template<class TImage>
  static bool SaveImages(std::vector<std::string> &filenames, const std::vector< typename itk::SmartPointer<TImage> > images);

  ///DICOM Related
  /*!
  \fn File::GetDICOMSeriesUIDs(const std::string directoryPath, bool recursive = false)
  \brief Returns the list of UIDs for a given directory.

  This is useful because directories may contain more than one dataset. Recursive allows for multiple datasets.

  Returns an empty list if an error was encountered.
  */
  static std::vector<std::string> GetDICOMSeriesUIDs(const std::string directoryPath, bool recursive = false);
  /*!
  \fn File::GetDICOMSeriesFilenames(const std::string directoryPath, const std::string seriesName)
  \brief Returns the filenames for a given UID/Series name for a given directory.

  This function simply returns only filenames and nothing is actually read as an image.

  Returns an empty list if an error was encountered.
  */
  static std::vector<std::string> GetDICOMSeriesFilenames(const std::string directoryPath, const std::string seriesName);
  /*!
  \fn File::GetDICOMData(const std::string directoryPath, typename itk::SmartPointer<TImage> &data, itk::MetaDataDictionary &dict, std::string &seriesName, std::string &caseID)
  \brief Opens a DICOM series and returns the image meta data read by ITK/GDCM, i.e. the DICOM tags from the header, as well as the image data.

  You can get the UIDs via the GetDICOMSeriesUIDs() member, which is read from the filenames. This does not read the whole image data and is faster for a quick query.
  */
  template<class TImage>
  static bool GetDICOMData(const std::string directoryPath, typename itk::SmartPointer<TImage> &data, itk::MetaDataDictionary &dict, std::string &seriesName, std::string &caseID);
  /*!
  \fn File::GetDICOMMetaData(const std::string directoryPath, itk::MetaDataDictionary &dict, std::string &seriesName, std::string &caseID = "")
  \brief Opens a DICOM series and returns the image meta data read by ITK/GDCM, i.e. the DICOM tags from the header.

  You can get the UIDs via the GetDICOMSeriesUIDs() member, which is read from the filenames. This does not read the whole image data and is faster for a quick query.
  */
  template<class TImage>
  static bool GetDICOMMetaData(const std::string directoryPath, itk::MetaDataDictionary &dict, std::string &seriesName, std::string &caseID);
  /*!
  \fn File::GetDICOMTags(const std::string directoryPath, std::vector< std::pair<std::string, std::string> > &tags, std::string &seriesName, std::string &caseID)
  \brief Opens a DICOM series and returns the DICOM tags read by ITK/GDCM, i.e. the DICOM tags from the header, as a list of string pairs.

  You can get the UIDs via the GetDICOMSeriesUIDs() member and the meta data via GetDICOMMetaData(), which are read from the filenames.
  */
  template<class TImage>
  static bool GetDICOMTags(const std::string directoryPath, std::vector< std::pair<std::string, std::string> > &tags, std::string &seriesName, std::string &caseID);
  template<class TImage>
  static bool GetDICOMTags(const itk::MetaDataDictionary dictionary, std::vector< std::pair<std::string, std::string> > &tags);
  /*!
  \fn File::OpenDICOMSeries(const std::string directoryPath, typename itk::SmartPointer<TImage> &data, std::string &seriesName, std::string &caseID)
  \brief Opens a DICOM series from the path given. Returns the image volume from the series read by ITK/GDCM.

  You can get the UIDs via the GetDICOMSeriesUIDs() member, which is read from the filenames. The DICOM tags are read and seriesName and caseID is replaced with DICOM tag values.

  Returns true if successful. Image is also NOT flipped for VTK.
  */
  template<class TImage>
  static bool OpenDICOMSeries(const std::string directoryPath, typename itk::SmartPointer<TImage> &data, std::string &seriesName, std::string &caseID);
  template<class TImage>
  static bool OpenDICOMSeriesAndTags(const std::string directoryPath, typename itk::SmartPointer<TImage> &data, std::vector< std::pair<std::string, std::string> > &tags, std::string &seriesName, std::string &caseID);
#endif

#ifndef ITK_ONLY
  #ifndef VTK_ONLY
  /*!
    \fn File::OpenImage(const std::string filename, vtkSmartPointer<vtkImageData> &data)
    \brief Opens an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI, INR etc. Overloaded for VTK input/output.

    ITK Reader is used and then transparently converted to a VTK image.
    Returns true if successful. Image is also flipped since the ITK reader orientation is different to VTK images.
  */
  template<class TImage>
  static bool OpenImage(const std::string filename, vtkSmartPointer<vtkImageData> &data);
  /*!
    \fn File::SaveImage(const std::string filename, vtkSmartPointer<vtkImageData> &data)
    \brief Saves an image file, which is any of the following: JPEG, PNG, DICOM, TIFF, NIFTI, INR etc. Overloaded for VTK input/output.

    ITK Writer is used after transparently converting from a VTK image.
    Returns true if successful. Image is also flipped since the ITK reader orientation is different to VTK images.
  */
  template<class TImage>
  static bool SaveImage(const std::string filename, vtkSmartPointer<vtkImageData> &data);
  #endif
#endif

#ifndef VTK_ONLY
  /**
    \brief Opens the image using the ITK file reader class. Returns NULL if failed and outputs the error to std error.

    This member is inline deliberately to avoid function call overheads.
  */
  template<class TImage>
  inline static itk::SmartPointer<TImage> ReadImageUsingITK(const std::string filename);
  /**
    \brief Saves the image using the ITK file writer class. Returns NULL if failed and outputs the error to std error.

    This member is inline deliberately to avoid function call overheads.
  */
  template<class TImage>
  inline static bool WriteImageUsingITK(const std::string filename, itk::SmartPointer<TImage> data);
#endif

  //Transforms
#ifndef VTK_ONLY
  /*!
      \fn File::OpenTransform(std::string filename)
      \brief Opens an ITK transform file and returns.
      \author Shekhar Chandra, 2013
  */
  template<class TType>
  static itk::SmartPointer< itk::Transform<TType> > OpenTransform(std::string filename);
#endif
#ifndef ITK_ONLY
  #ifndef VTK_ONLY
  /*!
      \fn File::OpenITKTransform(std::string filename)
      \brief Opens an ITK transform file and converts it to a VTK transform matrix.
      \author Jurgen Fripp et al., 2010. Modified by Shekhar Chandra, 2011
  */
  static vtkMatrix4x4* OpenITKTransform(std::string filename);

  /*!
      \fn File::OpenVTKTransform(std::string filename)
      \brief Opens a VTK transform file and converts it to a VTK transform matrix.
      \author Jurgen Fripp et al., 2010. Modified by Shekhar Chandra, 2011
  */
  static vtkMatrix4x4* OpenVTKTransform(std::string filename);

  /*!
      \fn milxFile::SaveITKTransform(std::string filename, const vtkMatrix4x4 * matrix)
      \brief Saves a VTK transform matrix in an ITK transform file.
      \author David Rivest-Henault, 2013
  */
  static void SaveITKTransform(std::string filename, const vtkMatrix4x4 * matrix);

  /*!
      \fn milxFile::SaveVTKTransform(std::string filename, const vtkMatrix4x4 * matrix)
      \brief Saves a VTK transform matrix in an ITK transform file. The recommended filename
      extension for the VTK transform file is '.trsf'. The file syntax is like, eg:
      (
      O8
      1.000000 0.000000  0.000000  22.500000
      0.000000 1.000000  0.000000  22.500000
      0.000000 0.000000  1.000000  22.500000
      0.000000 0.000000  0.000000  1.000000
      )
      Note: 08 is a magic number.
      \author David Rivest-Henault, 2013
  */
  static void SaveVTKTransform(std::string filename, const vtkMatrix4x4 * matrix);

  #endif
#endif

#ifndef ITK_ONLY
  //Surfaces/Models
  /*!
    \fn File::OpenModel(const std::string filename, vtkSmartPointer<vtkPolyData> &data)
    \brief Opens a model file, which can be a VTK XML, Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk), a Polygonal File (*.ply) or a Object file (*.obj).

    Returns true if successful. PolyData is allocated within this member, so pass a NULL pointer. Hence, the pass is by reference.
  */
  static bool OpenModel(const std::string filename, vtkSmartPointer<vtkPolyData> &data);
  /*!
    \fn File::OpenModelCollection(std::vector<std::string> filenames, vtkSmartPointer<vtkPolyDataCollection> &collection)
    \brief Opens model files, which can be a VTK XML, Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk), a Polygonal File (*.ply) or a Object file (*.obj) into a collection.

    Returns true if successful. PolyData and PolyDataCollection are allocated within this member, so pass a NULL pointer. Hence, the pass is by reference.
  */
  static bool OpenModelCollection(std::vector<std::string> filenames, vtkSmartPointer<vtkPolyDataCollection> &collection);
  /*!
    \fn File::SaveModel(const std::string filename, vtkSmartPointer<vtkPolyData> data, const bool binary = false)
    \brief Saves a model as a file, which can be a VTK XML or Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk) or a Polygonal File (*.ply).

    Returns true if successful.
    */
  static bool SaveModel(const std::string filename, vtkSmartPointer<vtkPolyData> data, const bool binary = false);
  /*!
    \fn File::SaveModelCollection(std::vector<std::string> filenames, vtkSmartPointer<vtkPolyDataCollection> collection, const bool binary = false)
    \brief Saves model files, which can be a VTK XML, Legacy VTK PolyData File (i.e. either a *.vtp or *.vtk), a Polygonal File (*.ply) or a Object file (*.obj) from a collection.

    Returns true if successful.
  */
  static bool SaveModelCollection(std::vector<std::string> filenames, vtkSmartPointer<vtkPolyDataCollection> collection, const bool binary = false);
  /*!
   \fn File::SaveTransform(const std::string filename, vtkSmartPointer<vtkTransform> data, const bool ITK = false)
   \brief Saves a transform file, which can be VTK .trsf or ITK from a collection.

   Returns true if successful.
   */
 static bool SaveTransform(const std::string filename, vtkSmartPointer<vtkTransform> data, const bool ITK = false);
  /*!
    \fn File::SaveTransformCollection(std::vector<std::string> filenames, vtkSmartPointer<vtkTransformCollection> collection, const bool ITK = false)
    \brief Saves transform files, which can be VTK .trsf or ITK from a collection.

    Returns true if successful.
  */
  static bool SaveTransformCollection(std::vector<std::string> filenames, vtkSmartPointer<vtkTransformCollection> collection, const bool ITK = false);

  //CSV files
  /*!
    \fn File::OpenDelimitedText(const std::string filename, vtkSmartPointer<vtkTable> &data, const std::string delimiter = ",")
    \brief Opens ASCII Delimited text files such as CSV files. Data is loaded into the table.

    Returns true if successful.
  */
  static bool OpenDelimitedText(const std::string filename, vtkSmartPointer<vtkTable> &data, const std::string delimiter = ",");
  /*!
    \fn File::SaveDelimitedText(const std::string filename, const vtkSmartPointer<vtkTable> data, const std::string delimiter = ",")
    \brief Saves table as an ASCII Delimited text files such as a CSV file.

    Returns true if successful.
  */
  static bool SaveDelimitedText(const std::string filename, const vtkSmartPointer<vtkTable> data, const std::string delimiter = ",");

  //Camera
  static bool SaveCamera(const std::string filename, const vtkSmartPointer<vtkCamera> camera);
  static bool LoadCamera(const std::string filename, vtkSmartPointer<vtkCamera> &camera);
#endif

  /*!
    \fn File::StringToLowerCase(std::string filename)
    \brief Returns a lower case version of the string.
  */
  static inline std::string StringToLowerCase(std::string filename)
  {
    std::transform(filename.begin(), filename.end(), filename.begin(), ::tolower);
    return filename;
  }
  /*!
    \fn File::StripFileExtension(const std::string filename)
    \brief Returns filename with last extension removed.

    "/foo/bar/baz.txt" --> "txt"
  */
  static inline std::string StripFileExtension(const std::string &filename)
  {
    return filename.substr(0, filename.find_last_of("."));
  }
  /*!
    \fn File::ExtractFilename(const std::string &filename, char delimiter = '/')
    \brief Returns filename (with extension) with paths removed.

    e.g. "/foo/bar/baz.txt" --> "baz.txt"
  */
  static inline std::string ExtractFilename(const std::string &filename, char delimiter = '/')
  {
    return filename.substr( filename.find_last_of( delimiter ) + 1 );
  }
  /*!
    \fn File::ExtractPath(const std::string &filename, char delimiter = '/')
    \brief Returns the path of the filename with the basename removed.

    "/foo/bar/baz.txt" --> "/foo/bar/"
  */
  static inline std::string ExtractPath(const std::string &filename, char delimiter = '/')
  {
    return filename.substr( 0, filename.find_last_of( delimiter ) + 1 );
  }
  /*!
    \fn File::GetBaseName(std::string filename)
    \brief Returns filename with paths and extension removed.

    e.g. "/foo/bar/baz.txt" --> "baz"
  */
  static inline std::string GetBaseName(const std::string &filename)
  {
    return StripFileExtension( ExtractFilename(filename) );
  }
  /*!
    \fn File::GetFileExtension(const std::string &filename)
    \brief Returns the file extension (in lower case) of the given filename

    This uses the last occurance of the period to get the right most substring.
  */
  static inline std::string GetFileExtension(const std::string &filename)
  {
    if (filename.find_last_of(".") != std::string::npos)
      return StringToLowerCase( filename.substr(filename.find_last_of(".")+1) );
    return "";
  }
  /*!
    \fn File::Exists(const std::string filename)
    \brief Returns true if the file already exists
  */
  static bool Exists(const std::string filename)
  {
    struct stat buffer;
    if (stat(filename.c_str(), &buffer) != -1)
        return true;
    return false;
  }

protected:

private:

};

#ifndef VTK_ONLY
template<class TImage>
bool File::OpenImage(const std::string filename, typename itk::SmartPointer<TImage> &data)
{
  if(!Exists(filename))
  {
    std::cerr << "File " << filename << " doesn't exist. Ignoring." << std::endl;
    return false;
  }

  data = ReadImageUsingITK<TImage>(filename);

  if(!data)
    return false;

  return true;
}

#if (ITK_VERSION_MAJOR > 3)
  template<class TImage, size_t TDim>
  //template<class TPixel, class TDim>
  bool File::OpenAsVectorImage(const std::string filename, typename itk::SmartPointer< itk::VectorImage< typename TImage::InternalPixelType, TDim> > &data)
  {
    if(!Exists(filename))
    {
      std::cerr << "File " << filename << " doesn't exist. Ignoring." << std::endl;
      return false;
    }

    ///Read the Hyper volume
    typename TImage::Pointer image = ReadImageUsingITK<TImage>(filename);

    if(!image)
      return false;

    typedef itk::Image< typename TImage::InternalPixelType, TDim > TImageSlice;
    typedef itk::ComposeImageFilter<TImageSlice> ImageToVectorImageFilterType;
    typename ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
    typename TImage::SizeType size = image->GetLargestPossibleRegion().GetSize();
    std::cout << "Image Size: " << size[0] << "x" << size[1] << "x" << size[2] << "x" << size[3] << std::endl;
    for(size_t j = 0; j < size[3]; j ++)
    {
      typename TImage::IndexType desiredStart; ///Pick the hyper plane
      desiredStart.Fill(0);
      desiredStart[3] = j;

      typename TImage::SizeType desiredSize;
      desiredSize[0] = size[0];
      desiredSize[1] = size[1];
      desiredSize[2] = size[2];
      desiredSize[3] = 0; //Needs to be zero to create slice

      typename TImage::RegionType desiredRegion(desiredStart, desiredSize);

      ///Extract hype plane, eg. a 3D image 'slice' of a 4D volume
      typedef itk::ExtractImageFilter<TImage, TImageSlice> FilterType;
      typename FilterType::Pointer filter = FilterType::New();
      filter->SetExtractionRegion(desiredRegion);
      filter->SetInput(image);
    #if ITK_VERSION_MAJOR >= 4
      filter->SetDirectionCollapseToIdentity(); // This is required.
    #endif
      filter->Update();

      ///add to vector image as a component
      imageToVectorImageFilter->SetInput(j, filter->GetOutput());
    }
    imageToVectorImageFilter->Update();

    data = imageToVectorImageFilter->GetOutput();

    if(!data)
      return false;

    return true;
  }
#endif

template<class TImage>
bool File::OpenImages(std::vector<std::string> &filenames, std::vector< typename itk::SmartPointer<TImage> > &images)
{
  if(filenames.empty())
    return false;

  std::vector<std::string> readFilenames;
  for(size_t j = 0; j < filenames.size(); j ++)
  {
    typename itk::SmartPointer<TImage> image;
    bool success = OpenImage(filenames[j], image);
    if(success && image)
    {
      readFilenames.push_back(filenames[j]);
      images.push_back(image);
    }
  }

  if(images.empty())
    return false;
  filenames = readFilenames;

  return true;
}

template<class TImage>
bool File::SaveImage(const std::string filename, typename itk::SmartPointer<TImage> data, itk::ImageIOBase *io)
{
  typedef itk::ImageFileWriter<TImage> ImageWriter;

  typename ImageWriter::Pointer writer = ImageWriter::New();
    writer->UseInputMetaDataDictionaryOn();
    writer->SetInput(data);
    writer->SetFileName(filename.c_str());
    if(io)
      writer->SetImageIO(io);
    writer->AddObserver(itk::ProgressEvent(), ProgressUpdates);
    try
    {
      writer->Update();
    }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "File Exception caught while writing!" << std::endl;
      std::cerr << err << std::endl;
      return false;
    }

  return true;
}

template<class TImage>
bool File::SaveImages(std::vector<std::string> &filenames, const std::vector< typename itk::SmartPointer<TImage> > images)
{
  if(filenames.empty())
    return false;

  std::vector<std::string> readFilenames;
  for(size_t j = 0; j < filenames.size(); j ++)
  {
    bool success = SaveImage<TImage>(filenames[j], images[j]);
    if(success && images[j])
      readFilenames.push_back(filenames[j]);
  }

  if(images.empty())
    return false;
  filenames = readFilenames;

  return true;
}

//DICOM related
template<class TImage>
bool File::GetDICOMData(const std::string directoryPath, typename itk::SmartPointer<TImage> &data, itk::MetaDataDictionary &dict, std::string &seriesName, std::string &caseID)
{
  typedef itk::ImageSeriesReader<TImage> ReaderType;
  typedef itk::GDCMImageIO ImageIOType;

  const std::vector<std::string> filenames = GetDICOMSeriesFilenames(directoryPath, seriesName);

  ImageIOType::Pointer gdcmIO = ImageIOType::New();
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetImageIO(gdcmIO);
  reader->SetFileNames(filenames);
  reader->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject &excp)
  {
    std::cerr << "File Exception caught while reading series!" << std::endl;
    std::cerr << excp << std::endl;
    return false;
  }

  ///Print info
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
  std::string series_type_id("0008|103e");
  itk::MetaDataDictionary & dic = gdcmIO->GetMetaDataDictionary();
  itk::MetaDataDictionary::ConstIterator series_type_itr = dic.Find(series_type_id);

  std::string caseId = "0010|0020";
  itk::MetaDataDictionary::ConstIterator case_itr = dic.Find(caseId);

  std::string echoNumber = "0018|0086";
  itk::MetaDataDictionary::ConstIterator echoNumber_itr = dic.Find(echoNumber);

  MetaDataStringType::ConstPointer entryValue = dynamic_cast<const MetaDataStringType *>(series_type_itr->second.GetPointer());
  if (entryValue)
  {
    seriesName = entryValue->GetMetaDataObjectValue();
    std::cout << "Series: " << entryValue->GetMetaDataObjectValue() << std::endl;
  }
  MetaDataStringType::ConstPointer entryValue2 = dynamic_cast<const MetaDataStringType *>(case_itr->second.GetPointer());
  if (entryValue2)
  {
    caseID = entryValue2->GetMetaDataObjectValue();
    std::cout << "Case: " << entryValue2->GetMetaDataObjectValue() << std::endl;
  }
  MetaDataStringType::ConstPointer entryValue3;
  if (dic.HasKey(echoNumber))
  {
    entryValue3 = dynamic_cast<const MetaDataStringType *>(echoNumber_itr->second.GetPointer());
    std::cout << "Echo Number: " << entryValue3->GetMetaDataObjectValue() << std::endl;
  }

  ///Return Meta Data as dictionary
  dict = gdcmIO->GetMetaDataDictionary();
  ///Return image data
  data = reader->GetOutput();

  return true;
}

template<class TImage>
bool File::GetDICOMMetaData(const std::string directoryPath, itk::MetaDataDictionary &dict, std::string &seriesName, std::string &caseID)
{
  typename itk::SmartPointer<TImage> data;

  if (!GetDICOMData<TImage>(directoryPath, data, dict, seriesName, caseID))
    return false;

  return true;
}

template<class TImage>
bool File::GetDICOMTags(const std::string directoryPath, std::vector< std::pair<std::string, std::string> > &tags, std::string &seriesName, std::string &caseID)
{
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
  typedef itk::MetaDataDictionary DictionaryType;
  DictionaryType dictionary;
  if (!GetDICOMMetaData<TImage>(directoryPath, dictionary, seriesName, caseID))
    return false;

  DictionaryType::ConstIterator itr = dictionary.Begin();
  DictionaryType::ConstIterator end = dictionary.End();

  while (itr != end)
  {
    itk::MetaDataObjectBase::Pointer entry = itr->second;
    MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>(entry.GetPointer());
    std::string entryTag = itr->first;
    if (entryvalue) //convert MetaDataObject pair to string pair
    {
      std::pair<std::string, std::string> entryPair = std::make_pair(entryTag, entryvalue->GetMetaDataObjectValue());
      tags.push_back(entryPair);
    }
    ++itr;
  }

  return true;
}

template<class TImage>
bool File::GetDICOMTags(const itk::MetaDataDictionary dictionary, std::vector< std::pair<std::string, std::string> > &tags)
{
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
  typedef itk::MetaDataDictionary DictionaryType;

  DictionaryType::ConstIterator itr = dictionary.Begin();
  DictionaryType::ConstIterator end = dictionary.End();

  while (itr != end)
  {
    itk::MetaDataObjectBase::Pointer entry = itr->second;
    MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>(entry.GetPointer());
    std::string entryTag = itr->first;
    if (entryvalue) //convert MetaDataObject pair to string pair
    {
      std::pair<std::string, std::string> entryPair = std::make_pair(entryTag, entryvalue->GetMetaDataObjectValue());
      tags.push_back(entryPair);
    }
    ++itr;
  }

  return true;
}

template<class TImage>
bool File::OpenDICOMSeries(const std::string directoryPath, typename itk::SmartPointer<TImage> &data, std::string &seriesName, std::string &caseID)
{
  typedef itk::MetaDataDictionary DictionaryType;
  DictionaryType dictionary;

  if (!GetDICOMData<TImage>(directoryPath, data, dictionary, seriesName, caseID))
    return false;

  return true;
}

template<class TImage>
bool File::OpenDICOMSeriesAndTags(const std::string directoryPath, typename itk::SmartPointer<TImage> &data, std::vector< std::pair<std::string, std::string> > &tags, std::string &seriesName, std::string &caseID)
{
  typedef itk::MetaDataDictionary DictionaryType;
  DictionaryType dictionary;

  if (!GetDICOMData<TImage>(directoryPath, data, dictionary, seriesName, caseID))
    return false;

  if (!milx::File::GetDICOMTags<TImage>(dictionary, tags))
    return false;

  return true;
}
#endif

#ifndef ITK_ONLY
  #ifndef VTK_ONLY
  template<class TImage>
  bool File::OpenImage(const std::string filename, vtkSmartPointer<vtkImageData> &data)
  {
    std::string extension = GetFileExtension(filename);
    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();

    if(!Exists(filename))
      return false;

    if(extension == "vti") //VTK XML Image format
    {
      vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();

      if(reader->CanReadFile(filename.c_str()))
      {
          reader->SetFileName(filename.c_str());
          reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
          //~ linkProgressEventOf(reader);
          reader->Update();

          if(!errorObserver->ReportsFailure())
          {
              data = reader->GetOutput();
          }
          else
          {
              std::cerr << "VTK XML Image Reader Encountered the following error." << std::endl;
              std::cerr << errorObserver->GetMessage() << std::endl;
              return false;
          }
      }
      else
      {
          std::cerr << "Could not locate VTK XML Image (VTI) file!" << std::endl;
          return false;
      }
    }
    else
    {
      typename itk::SmartPointer<TImage> tmpImage;

      if( !OpenImage(filename, tmpImage) )
        return false;

      typedef itk::ImageToVTKImageFilter<TImage> ConnectorType;
      ///Export to VTK
      typename ConnectorType::Pointer connector = ConnectorType::New(); //!< ITK to VTK image
        connector->SetInput(tmpImage);
        connector->Update();

      ///Flip - ITK and VTK do not have the same orientation
      vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
      #if VTK_MAJOR_VERSION <= 5
        imageReorient->SetInput(connector->GetOutput());
      #else
        imageReorient->SetInputData(connector->GetOutput());
      #endif
        imageReorient->SetFilteredAxis(1);
        imageReorient->FlipAboutOriginOn();
        imageReorient->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        imageReorient->Update();

      if(errorObserver->ReportsFailure())
      {
        std::cerr << "Reader Encountered the following error." << endl;
        std::cerr << errorObserver->GetMessage() << endl;
        return false;
      }

      data = imageReorient->GetOutput();
    }

    return true;
  }

  template<class TImage>
  bool File::SaveImage(const std::string filename, vtkSmartPointer<vtkImageData> &data)
  {
    std::string extension = GetFileExtension(filename);

    ///Flip - ITK and VTK do not have the same orientation
    vtkSmartPointer<vtkErrorWarning> errorObserver = vtkSmartPointer<vtkErrorWarning>::New();
    vtkSmartPointer<vtkImageFlip> imageReorient = vtkSmartPointer<vtkImageFlip>::New();
    #if VTK_MAJOR_VERSION <= 5
      imageReorient->SetInput(data);
    #else
      imageReorient->SetInputData(data);
    #endif
      imageReorient->SetFilteredAxis(1);
      imageReorient->FlipAboutOriginOn();
      imageReorient->AddObserver(vtkCommand::ErrorEvent, errorObserver);
      imageReorient->Update();

    if(errorObserver->ReportsFailure())
    {
      std::cerr << "Writer Encountered the following error." << std::endl;
      std::cerr << errorObserver->GetMessage() << std::endl;
      return false;
    }

    if(extension == "vti") //VTK XML Image format
    {
      vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
      #if VTK_MAJOR_VERSION <= 5
        writer->SetInput(data);
      #else
        writer->SetInputData(data);
      #endif
        writer->SetFileName(filename.c_str());
        writer->SetDataModeToBinary();
        writer->AddObserver(vtkCommand::ErrorEvent, errorObserver);
        //~ linkProgressEventOf(writer);
        writer->Write();

      if(errorObserver->ReportsFailure())
      {
          std::cerr << "VTK XML Image Writer Encountered the following error." << std::endl;
          std::cerr << errorObserver->GetMessage() << std::endl;
      }
    }
    else
    {
      typedef itk::VTKImageToImageFilter<TImage> ConnectorType;
      ///Export to VTK
      typename ConnectorType::Pointer connector = ConnectorType::New(); //!< ITK to VTK image
        connector->SetInput(imageReorient->GetOutput());
        connector->Update();

      if( !SaveImage<TImage>(filename, connector->GetOutput()) )
        return false;
    }

    return true;
  }
  #endif
#endif

#ifndef VTK_ONLY
template<class TImage>
itk::SmartPointer<TImage> File::ReadImageUsingITK(const std::string filename)
{
  typedef itk::ImageFileReader<TImage, itk::DefaultConvertPixelTraits<typename TImage::InternalPixelType> > ImageReader; //InternalPixelType != PixelType for vector images
  typename ImageReader::Pointer reader = ImageReader::New();
  reader->SetFileName(filename.c_str());
  reader->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    reader->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "Reader Encountered the following error." << std::endl;
    std::cerr << err << std::endl;
    return NULL;
  }

  return reader->GetOutput();
}

template<class TImage>
bool File::WriteImageUsingITK(const std::string filename, itk::SmartPointer<TImage> data)
{
  typedef itk::ImageFileWriter<TImage> ImageWriter;
  typename ImageWriter::Pointer writer = ImageWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInput(data);
  writer->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "Writer Encountered the following error." << std::endl;
    std::cerr << err << std::endl;
    return false;
  }

  return true;
}

template<class TType>
itk::SmartPointer< itk::Transform<TType> > File::OpenTransform(std::string filename)
{
  typedef itk::TransformFileReader        TransformReaderType;
  typedef TransformReaderType::TransformListType * TransformListType;

  // Register default transforms
  itk::TransformFactoryBase::RegisterDefaultTransforms();

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
    return NULL;
  }

  return static_cast< itk::Transform<TType> *>( (*(affineReader->GetTransformList()->begin())).GetPointer() );
}
#endif

} //end namespace milx

#endif //__MILXFILE_H
