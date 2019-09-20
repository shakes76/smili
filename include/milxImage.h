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
#ifndef __MILXIMAGE_H
#define __MILXIMAGE_H
//ITK 4.13
#include <itkSpatialOrientationAdapter.h>
//ITK
#include <itkVectorImage.h>
#include <itkImageDuplicator.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkAdaptiveHistogramEqualizationImageFilter.h>
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
  #include <itkLabelImageToLabelMapFilter.h>
  #include <itkLabelMapToLabelImageFilter.h>
  #include <itkAutoCropLabelMapFilter.h>
  #include <itkBinaryContourImageFilter.h>
  #include <itkLabelContourImageFilter.h>
  #include <itkLabelOverlayImageFilter.h>
  #include <itkMergeLabelMapFilter.h>
#endif // (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
#include <itkImportImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkNormalizeImageFilter.h>
#include <itkInvertIntensityImageFilter.h>
#include <itkChangeInformationImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>
#include <itkCheckerBoardImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkApproximateSignedDistanceMapImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkOtsuMultipleThresholdsImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkFlipImageFilter.h>
#include <itkConstantPadImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkBilateralImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkSobelEdgeDetectionImageFilter.h>
#include <itkMedianImageFilter.h>
#include <itkLaplacianRecursiveGaussianImageFilter.h>
#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkExtractImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkAffineTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkCenteredEuler3DTransform.h>
#include <itkEuler3DTransform.h>
#include <itkTransformFileWriter.h>
#include <itkShrinkImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkJoinSeriesImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
//ITK4 Compatibility
#if (ITK_VERSION_MAJOR > 3)
  #include <itkv3Rigid3DTransform.h> //ITK4 (4.3.1) version has New() member issues
  #include <itkFFTConvolutionImageFilter.h>
  #include <itkVectorMagnitudeImageFilter.h> //vector images
#else
  #include <itkRigid3DTransform.h>
  #if (ITK_REVIEW || ITK_VERSION_MAJOR > 3) //Review only members
    #include <itkConvolutionImageFilter.h>
  #endif
  #include <itkGradientToMagnitudeImageFilter.h> //vector images
#endif
//Exporter/Importer
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#ifndef ITK_ONLY
  //VTK
  #include <vtkMatrix4x4.h>
  #include <vtkTransform.h>
  #include <vtkImageReslice.h>
#endif
//ITK Extensions
#include "itkVectorShiftScaleImageFilter.h"
//SMILI
#include "milxGlobal.h"

const double CoordinateTolerance = 1e-3;
const double DirectionTolerance = 1e-3;

namespace milx
{

/**
  \class ImageBase
  \brief Base class for images that contains non-templated features so to allow dynamic template instantiation.

  Please use the milx::Image class instead unless you require something like below:
  \code

  \endcode
*/
class SMILI_EXPORT ImageBase
{
public:
  /*!
  	\fn Image::ImageBase()
  	\brief Standard constructor
  */
  ImageBase();
  /*!
  	\fn Image::~ImageBase()
  	\brief Standard Destructor
  */
  virtual ~ImageBase() {};
};

/**
  \class Image
  \brief Represents an image (i.e. an regular rectangular array with scalar values) and their common operations using the Insight Toolkit (ITK).

  Due to the templated nature of itk::Image objects, most of the members in this class are static and the current state is not tracked.
  This class is used extensively throughout the milxQtImage class.
  Use the milx::File object to load and save images.
  See the usage examples given below:

  Converting
  \code
  typedef unsigned char InputPixelType;
  typedef itk::Image<InputPixelType, 3> LabelImageType;
  LabelImageType::SpacingType labelSpacing = m_LabelledImages[index]->GetSpacing();

  milx::Model model1(surface1);
  vtkSmartPointer<vtkImageData> voxelisedModel1 = model1.Voxelise(255, labelSpacing.GetDataPointer(), bigBounds);
  LabelImageType::Pointer image1 = milx::Image<LabelImageType>::ConvertVTKImageToITKImage(voxelisedModel1);
  \endcode

  Distance maps
  \code
  typedef float FloatPixelType;
  typedef itk::Image<FloatPixelType, 3> FloatImageType;
  const bool binary = false, signedDistance = true, insideDistance = false, squaredDistance = false;

  FloatImageType::Pointer distanceImage1 = milx::Image<LabelImageType>::DistanceMap<FloatImageType>(image1, binary, signedDistance, insideDistance, squaredDistance);
  \endcode

  Padding an image
  \code
  LabelImageType::RegionType region = image->GetLargestPossibleRegion();
  LabelImageType::SizeType imageSize = region.GetSize();
  itk::SmartPointer<LabelImageType> paddedImage = milx::Image<LabelImageType>::PadImageByConstant(image, imageSize[0]/2, imageSize[1]/2, imageSize[2]/2, 0);
  \endcode
*/
template<class TImage>
class SMILI_EXPORT Image : public ImageBase
{
public:
  /*!
  	\fn Image::Image()
  	\brief Standard constructor
  */
  Image();
  /*!
  	\fn Image::Image(itk::SmartPointer<TImage> image)
  	\brief Constructor that copies the input model.
  */
//  Image(itk::SmartPointer<TImage> image);
  /*!
  	\fn Image::~Image()
  	\brief Standard Destructor
  */
  virtual ~Image() {};

  /*!
  	\fn Image::SetInput(itk::SmartPointer<TImage> image)
  	\brief Assigns the input image to class
  */
//  void SetInput(itk::SmartPointer<TImage> image);
  /*!
  	\fn Image::Result()
  	\brief Returns the current image, i.e. the result of the latest operation.

  Could be NULL, the user must check.
  */
  inline itk::SmartPointer<TImage> Result()
  {
    return CurrentImage;
  }
  /*!
  	\fn Image::PreviousResult()
  	\brief Returns the previous image, i.e. the result of the penultimate operation.

  Could be NULL, the user must check.
  */
  inline itk::SmartPointer<TImage> PreviousResult()
  {
    return PreviousImage;
  }
  /*!
  	\fn Image::GetOutput()
  	\brief Returns the current image, i.e. the result of the latest operation ITK/VTK style.
  */
  inline itk::SmartPointer<TImage> GetOutput()
  {
    return Result();
  }

  /**
  * \defgroup Atomic Atomic Operations Members
  * \brief Members for atomic (one image at a time) operations
  */
  //@{
  //Conversions
  /*!
    \fn Image::ConvertITKImageToVTKImage(itk::SmartPointer<TImage> img)
    \brief Converts a ITK image object to an VTK image object. You MUST DeepCopy the result as the ITK smartpointer is not aware of the VTK smartpointer.
  */
  static vtkSmartPointer<vtkImageData> ConvertITKImageToVTKImage(itk::SmartPointer<TImage> img);
  /*!
    \fn Image::ConvertITKVectorImageToVTKImage(itk::SmartPointer<TImage> img)
    \brief Converts a ITK Vector image object to an VTK image object. You MUST DeepCopy the result as the ITK smartpointer is not aware of the VTK smartpointer.
  */
  static vtkSmartPointer<vtkImageData> ConvertITKVectorImageToVTKImage(itk::SmartPointer<TImage> img);

  template<typename TPrecision>
  static itk::SmartPointer<TImage> ApplyOrientationToITKImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> refImage, const bool labelledImage, const bool flipY = true, const bool ignoreDirection = false);
#ifndef ITK_ONLY //Requires VTK
  /*!
    \fn Image:: ApplyOrientationToVTKImage(vtkSmartPointer<vtkImageData> img, itk::SmartPointer<TImage> refImage, vtkSmartPointer<vtkMatrix4x4> &transformMatrix, const bool labelledImage, const bool flipY = true)
    \brief Applies orientation/direction and origin to a VTK image from a reference image.

    Flip y-axis to comply with VTK coordinate system? The transformMatrix, where the extraneous transforms (from the image orientation, such as flipping) will be stored, need not be pre-allocated.

    \code
    vtkSmartPointer<vtkImageData> newImageData = vtkSmartPointer<vtkImageData>::New();
    newImageData->DeepCopy( milx::Image<floatImageType>::ConvertITKImageToVTKImage(imageFloat) );
    imageData = milx::Image<floatImageType>::ApplyOrientationToVTKImage(newImageData, imageFloat, transformMatrix, true, true);
    \endcode
  */
  static vtkSmartPointer<vtkImageData> ApplyOrientationToVTKImage(vtkSmartPointer<vtkImageData> img, itk::SmartPointer<TImage> refImage, vtkSmartPointer<vtkMatrix4x4> &transformMatrix, const bool labelledImage, const bool flipY = true);
#endif
  /*!
    \fn Image::ConvertVTKImageToITKImage(vtkSmartPointer<vtkImageData> img)
    \brief Converts a VTK image object to an ITK image object.
  */
  static itk::SmartPointer<TImage> ConvertVTKImageToITKImage(vtkSmartPointer<vtkImageData> img);
  /*!
    \fn Image::CastImage(itk::SmartPointer<TImage> img)
    \brief Casts an image from one type to another.
  */
  template<typename TOutImage>
  static itk::SmartPointer<TOutImage> CastImage(itk::SmartPointer<TImage> img);
  /*!
    \fn Image::DeepCopy(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToCopyFrom)
    \brief Copy the contents of an image to another.
  */
  //static void DeepCopy(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToCopyFrom);
  /*!
    \fn Image::DuplicateImage(itk::SmartPointer<TImage> img)
    \brief Duplicates the image into a new image.
  */
  static itk::SmartPointer<TImage> DuplicateImage(itk::SmartPointer<TImage> img);
  /*!
    \fn Image::DuplicateImage(const TImage *img)
    \brief Duplicates the image into a new image.
  */
  static itk::SmartPointer<TImage> DuplicateImage(const TImage *img);
  /*!
    \fn Image::ImportVectorToImage(vnl_vector<TVector> &vec, typename TImage::SizeType size, itk::SmartPointer<TImage> image = NULL)
    \brief Imports a VNL vector to an ITK image object.

    Assumes the vector is not empty and that if the image is provided, it is setup correctly (spacing, origin etc.)
  */
  template<typename TVector>
  static itk::SmartPointer<TImage> ImportVectorToImage(vnl_vector<TVector> &vec, typename TImage::SizeType size, itk::SmartPointer<TImage> image = NULL);
  /*!
    \fn Image::ImportMatrixToImage(vnl_matrix<TMatrix> &matrix, itk::SmartPointer<TImage> image = NULL)
    \brief Imports a VNL matrix to an ITK image object.

    Assumes the vector is not empty and that if the image is provided, it is setup correctly (spacing, origin etc.).
  */
  template<typename TMatrix>
  static itk::SmartPointer<TImage> ImportMatrixToImage(vnl_matrix<TMatrix> &matrix, itk::SmartPointer<TImage> image = NULL);

  /*!
    \fn Image::BlankImage(const TImage::PixelType value, const TImage::SizeType imgSize)
    \brief Creates an empty image filled with value given.
  */
  static itk::SmartPointer<TImage> BlankImage(typename TImage::PixelType value, typename TImage::SizeType imgSize);

  //Geometry
  /*!
    \fn Image::ExtractSubImage(itk::SmartPointer<TImage> img, typename TImage::RegionType imgRegion)
    \brief Extract a sub image (such as a slice) from the current image given the region.
  */
  template<typename TSubImage>
  static itk::SmartPointer<TSubImage> ExtractSubImage(itk::SmartPointer<TImage> img, typename TImage::RegionType imgRegion);
  /*!
    \fn Image::ExtractSlice(itk::SmartPointer<TImage> img, int *extent)
    \brief Extract a sub image (such as a slice) from the current image given the extent (typically obtained from a VTK image actor etc.).

    Note that extent must be of length 2*dimension. No checks are done.
  */
  template<typename TImageSlice>
  static itk::SmartPointer<TImageSlice> ExtractSlice(itk::SmartPointer<TImage> img, int *extent);
  /*!
    \fn Image::ExtractComponent(itk::SmartPointer<TImage> img, int component)
    \brief Extract a component from the current vector image given the component index.

    No checks are done.
  */
  template<typename TImageComponent>
  static itk::SmartPointer<TImageComponent> ExtractComponent(itk::SmartPointer<TImage> img, int component);
  /*!
    \fn Image::ResizeImage(itk::SmartPointer<TImage> img, typename TImage::SizeType imgSize)
    \brief Resizes current image using current spacing.
  */
  static itk::SmartPointer<TImage> ResizeImage(itk::SmartPointer<TImage> img, typename TImage::SizeType imgSize);
  /*!
    \fn Image::ResizeImage(itk::SmartPointer<TImage> img, typename TImage::SizeType imgSize)
    \brief Resizes current image using given spacing.
  */
  static itk::SmartPointer<TImage> ResizeImage(itk::SmartPointer<TImage> img, typename TImage::SizeType imgSize, typename TImage::SpacingType outputSpacing);
  /*!
    \fn Image::ResizeImage(itk::SmartPointer<TImage> img, typename TImage::SizeType imgSize, typename TImage::SpacingType outputSpacing, typename TImage::PointType outputOrigin, typename TImage::DirectionType outputDirection)
    \brief Resizes current image using given spacing, origin and direction.
  */
  static itk::SmartPointer<TImage> ResizeImage(itk::SmartPointer<TImage> img, typename TImage::SizeType imgSize, typename TImage::SpacingType outputSpacing, typename TImage::PointType outputOrigin, typename TImage::DirectionType outputDirection);
  /*!
    \fn Image::SubsampleImage(itk::SmartPointer<TImage> img, const size_t *factors)
    \brief Subsamples or shrinks the current image using the factors provided.

    Using a factor of 2 reduces the size of the dimension in the image by half.
    factors must match image dimension, this is why SizeType is used.
  */
  static itk::SmartPointer<TImage> SubsampleImage(itk::SmartPointer<TImage> img, typename TImage::SizeType factors);

  //Transform
  /*!
    \fn Image::TransformImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TOutImage> refImg, itk::SmartPointer<TTransform> transf, const bool inverse, const int interp = 1)
    \brief Transform the image into a new reference image space given the transform.

    Supports Affine, Versor Rigid 3D, Rigid 3D, Centered Euler 3D and Euler 3D transforms
    Interpolation variable: 0 - NN, 1 - Linear, 2 - BSpline
  */
  template<typename TOutImage, typename TTransform, typename TPrecision>
  static itk::SmartPointer<TOutImage> TransformImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TOutImage> refImg, itk::SmartPointer<TTransform> transf, const bool inverse, const int interp = 1);
  template<typename TOutImage, typename TTransform, typename TPrecision>
  static itk::SmartPointer<TOutImage> TransformImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TTransform> transf, const bool inverse, const int interp = 1);
  /*!
    \fn Image::ResampleImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TOutImage> refImg, const bool linearInterp = false)
    \brief Resample the image into a new reference image space given.
  */
  template<typename TOutImage>
  static itk::SmartPointer<TOutImage> ResampleImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TOutImage> refImg, const bool linearInterp = false);
  /*!
    \fn Image::ResampleLabel(itk::SmartPointer<TImage> img, itk::SmartPointer<TOutImage> refImg)
    \brief Resample the a labelled image into a new reference image space given ensuring no interpolator artefacts.
  */
  template<typename TOutImage>
  static itk::SmartPointer<TOutImage> ResampleLabel(itk::SmartPointer<TImage> img, itk::SmartPointer<TOutImage> refImg);

  //Arithmetic
  /*!
    \fn Image::AddImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2)
    \brief Adds image 1 to image 2, returning the result.
  */
  static itk::SmartPointer<TImage> AddImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2);
  /*!
    \fn Image::SubtractImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2)
    \brief Subtracts image 2 from image 1, returning the result.
  */
  static itk::SmartPointer<TImage> SubtractImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2);
  /*!
    \fn Image::DifferenceImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2)
    \brief Same as SubtractImages().
  */
  static inline itk::SmartPointer<TImage> DifferenceImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2)
  {
    return SubtractImages(img1, img2);
  }
  /*!
  \fn Image::MultiplyImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2)
  \brief Multiplies (element-wise) image 2 from image 1, returning the result.
  */
  static itk::SmartPointer<TImage> MultiplyImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2);
  /*!
    \fn Image::ScaleImage(itk::SmartPointer<TImage> img, float scaling)
    \brief Scales the image intensities by scaling factor and returns the result.
  */
  template<class TOutImage>
  static itk::SmartPointer<TOutImage> ScaleImage(itk::SmartPointer<TImage> img, float scaling);
  /*!
    \fn Image::ScaleVectorImage(itk::SmartPointer<TImage> img, float scaling, int numberOfComponents)
    \brief Scales each component of the vector image intensities by scaling factor and returns the result.
  */
  static itk::SmartPointer<TImage> ScaleVectorImage(itk::SmartPointer<TImage> img, float scaling, int numberOfComponents);
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
  /*!
    \fn Image::ConvolveImages(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> kernelImg)
    \brief Convolves two images together and returns the result.

    If ITK 4 is found then the convolution theorem is used, i.e. convolution is computed in \f$ N\logN \f$, where \f$ N^n \f$ is the size of the image and \f$ n \f$ is the dimension.
  */
  static itk::SmartPointer<TImage> ConvolveImages(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> kernelImg);
#endif // (ITK_REVIEW || ITK_VERSION_MAJOR > 3)

  //Info
  /*!
  \fn Image::Information(itk::SmartPointer<TImage> img)
  \brief Prints the information of the image to standard output.
  */
  static void Information(itk::SmartPointer<TImage> img);
  /*!
  \fn Image::ImageMaximum(itk::SmartPointer<TImage> img)
  \brief Returns the maximum pixel value of an image.
  */
  static double ImageMaximum(itk::SmartPointer<TImage> img);
  /*!
  \fn Image::ImageMinimum(itk::SmartPointer<TImage> img)
  \brief Returns the minimum pixel value of an image.
  */
  static double ImageMinimum(itk::SmartPointer<TImage> img);
  /*!
  \fn Image::ImageOrientation(itk::SmartPointer<TImage> img)
  \brief Returns the orientation flag of an image.
  */
  static std::string ImageOrientation(itk::SmartPointer<TImage> img);

  //Filters
  /*!
  	\fn Image::CheckerBoard(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToCheckerBoard, const int numberOfSquares = 0)
  	\brief Generates a checker board pattern where each box alternates between two images. Ideal for comparison among images.
  */
  static itk::SmartPointer<TImage> CheckerBoard(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToCheckerBoard, const int numberOfSquares = 0);
  /*!
  	\fn Image::RescaleIntensities(itk::SmartPointer<TImage> img, float min, float max)
  	\brief Generates an image with the intensities rescaled as directed.
  */
  static itk::SmartPointer<TImage> RescaleIntensity(itk::SmartPointer<TImage> img, float minValue, float maxValue);
  /*!
  	\fn Image::InvertIntensity(itk::SmartPointer<TImage> img, float maxValue)
  	\brief Generates an image with the intensities reversed based on the max pixel value.
  */
  static itk::SmartPointer<TImage> InvertIntensity(itk::SmartPointer<TImage> img, float maxValue);
  /*!
  	\fn Image::HistogramEqualisation(itk::SmartPointer<TImage> img, float alpha, float beta)
  	\brief Generates an image with the intensities after histogram equalisation. Defaults to classic histogram equalisation.

  	Use parameters to change equalisation type/style. See reference
  	"Adaptive Image Contrast Enhancement using Generalizations of Histogram Equalization."
  	J.Alex Stark. IEEE Transactions on Image Processing, May 2000.
  */
  static itk::SmartPointer<TImage> HistogramEqualisation(itk::SmartPointer<TImage> img, float alpha = 0.3, float beta = 0.3, float radius = 5);
  /*!
  	\fn Image::GradientMagnitude(itk::SmartPointer<TImage> img)
  	\brief Generates an image with the gradient magnitude of the given image.
  */
  static itk::SmartPointer<TImage> GradientMagnitude(itk::SmartPointer<TImage> img);
  /*!
  	\fn Image::SobelEdges(itk::SmartPointer<TImage> img)
  	\brief Generates an image with Sobel edges, i.e. uses the Sobel edge detection on the input image.
  */
  static itk::SmartPointer<TImage> SobelEdges(itk::SmartPointer<TImage> img);
  /*!
  	\fn Image::Laplacian(itk::SmartPointer<TImage> img)
  	\brief Generates an Laplacian image from the input image.
  */
  static itk::SmartPointer<TImage> Laplacian(itk::SmartPointer<TImage> img);
  /*!
  	\fn Image::CannyEdges(itk::SmartPointer<TImage> img, float variance, float lowerThreshold, float upperThreshold)
  	\brief Generates an Canny edge image from the input image.
  */
  static itk::SmartPointer<TImage> CannyEdges(itk::SmartPointer<TImage> img, float variance, float lowerThreshold, float upperThreshold);
  /*!
  	\fn Image::Normalization(itk::SmartPointer<TImage> img)
  	\brief Generates an image which is statistically normalised so having pixel values between -1 and 1.
  */
  static itk::SmartPointer< itk::Image<float, TImage::ImageDimension> > Normalization(itk::SmartPointer<TImage> img);

  /*!
  	\fn Image::MatchInformation(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToMatch, bool originOnly = false)
  	\brief Changes the image info to match that of provided image. By default, only the spacing, region and origin are changed.

  	Set the originOnly argument true to only match the origin.
  */
  static itk::SmartPointer<TImage> MatchInformation(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToMatch, bool originOnly = false);
  /*!
    \fn Image::CopyInformation(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToMatch, bool originOnly = false)
    \brief Changes the image info to match that of provided image. By default, only the spacing, region and origin are changed.

    Set the originOnly argument true to only match the origin. Function is alias for MatchInformation().
  */
  static itk::SmartPointer<TImage> CopyInformation(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToMatch, bool originOnly = false)
  {
    return MatchInformation(img, imgToMatch, originOnly);
  }

  /*!
  	\fn Image::MatchHistogram(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToMatch, int bins = 128)
  	\brief Changes the image gray levels to match histogram of image provided.
  */
  static itk::SmartPointer<TImage> MatchHistogram(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToMatch, const int bins = 128);
  /*!
  	\fn Image::DistanceMap(itk::SmartPointer<TImage> img, const bool binaryImage = false, const bool signedDistances = false, const bool computeInsideObject = false, const bool squaredDistance = false)
  	\brief Generates a distance map image, where the distances are from the object boundary to the outside or vice versa.

  	The inside is considered as having negative distances. Outside is treated as having positive distances.
  	Use the signed argument to get signed distances using the Maurer method (TPAMI, 2003). Using the computeInside argument to compute distance within the object rather than outside.
  	If binary image is set to true, integer only computation is done if possible.
  */
  template<typename TOutImage>
  static itk::SmartPointer<TOutImage> DistanceMap(itk::SmartPointer<TImage> img, const bool binaryImage = false, const bool signedDistances = false, const bool computeInsideObject = false, const bool squaredDistance = false);
  /*!
  	\fn Image::FlipImage(itk::SmartPointer<TImage> img, bool xAxis = false, bool yAxis = true, bool zAxis = false, bool aboutOrigin = true)
  	\brief Generates an image with axes flipped. Currently flips the y-axis as default, but can do other axes also in combination.
  */
  static itk::SmartPointer<TImage> FlipImage(itk::SmartPointer<TImage> img, bool xAxis = false, bool yAxis = true, bool zAxis = false, bool aboutOrigin = true);
  /*!
  	\fn Image::PadImageByConstant(itk::SmartPointer<TImage> img, size_t xAxis, size_t yAxis, size_t zAxis, typename TImage::PixelType value)
  	\brief Generates an image padded by extending each axes by size given in both directions. The created areas have the value given by value.
  */
  static itk::SmartPointer<TImage> PadImageByConstant(itk::SmartPointer<TImage> img, size_t xAxis, size_t yAxis, size_t zAxis, typename TImage::PixelType value);
  /*!
  	\fn Image::AnisotropicDiffusion(itk::SmartPointer<TImage> img, const int iterations, const float timestep)
  	\brief Generates the gradient anisotropic diffusion (smoothing) of an image using the number of iterations.
  */
  template<typename TOutImage>
  static itk::SmartPointer<TOutImage> AnisotropicDiffusion(itk::SmartPointer<TImage> img, const int iterations, const float timestep);
  /*!
    \fn Image::Bilateral(itk::SmartPointer<TImage> img, const float sigmaRange, const float sigmaSpatial)
    \brief Generates the (non-linear) Bilateral smoothing of an image using the sigmas provided.
  */
  static itk::SmartPointer<TImage> Bilateral(itk::SmartPointer<TImage> img, const float sigmaRange, const float sigmaSpatial);
  /*!
    \fn Image::GaussianSmooth(itk::SmartPointer<TImage> img, const float variance)
    \brief Generates the Gaussian smoothing (by convolution) of an image using the variance provided.
  */
  static itk::SmartPointer<TImage> GaussianSmooth(itk::SmartPointer<TImage> img, const float variance);
  /*!
  	\fn Image::Median(itk::SmartPointer<TImage> img, const int radius)
    \brief Generates the (non-linear) median filtering (smoothing) of an image of given neighbourhood radius.
  */
  static itk::SmartPointer<TImage> Median(itk::SmartPointer<TImage> img, const int radius);

  //Labelling
  /*!
  	\fn Image::MaskImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> maskImg)
  	\brief Mask an image by given binary image.
  */
  template<typename TMaskImage>
  static itk::SmartPointer<TImage> MaskImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TMaskImage> maskImg);
  /*!
  	\fn Image::RelabelImage(itk::SmartPointer<TImage> labelledImg)
  	\brief Returns a labelled image with labelled values relabelled consecutively based on connectivity.
  */
  static itk::SmartPointer<TImage> RelabelImage(itk::SmartPointer<TImage> labelledImg);
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
  /*!
  	\fn Image::BinaryContour(itk::SmartPointer<TImage> img, const float foregroundValue, const float backgroundValue)
  	\brief Generates a contour from a binary image and returns it.
  */
  static itk::SmartPointer<TImage> BinaryContour(itk::SmartPointer<TImage> img, const float foregroundValue, const float backgroundValue);
  /*!
  	\fn Image::LabelContour(itk::SmartPointer<TImage> img, const bool fullyConnected, const float backgroundValue)
  	\brief Generates a contour for each label in the labelled image and returns it.
  */
  static itk::SmartPointer<TImage> LabelContour(itk::SmartPointer<TImage> img, const bool fullyConnected, const float backgroundValue);
  /*!
  	\fn Image::MaskAndCropImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> maskImg, TImage::SizeType::SizeValueType pixelPadding = 1)
  	\brief Mask and crop an image by given binary image using its bounds padded by given number of pixels.
  */
  template<typename TMaskImage>
  static itk::SmartPointer<TImage> MaskAndCropImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TMaskImage> maskImg, const size_t pixelPadding = 1);
  /*!
  	\fn Image::LabelValues(itk::SmartPointer<TImage> img)
  	\brief Returns a list of labelled values within the labelled image.
  */
  static std::vector<unsigned char> LabelValues(itk::SmartPointer<TImage> labelledImg);
  /*!
  	\fn Image::Overlay(itk::SmartPointer<TImage> img, itk::SmartPointer<TLabelImage> overlayImg)
  	\brief Generates an overlay of the images provided. Returns an RGB image of the result.
  */
  template<typename TLabelImage>
  static itk::SmartPointer< itk::Image<itk::RGBPixel<unsigned char>, TImage::ImageDimension> > Overlay(itk::SmartPointer<TImage> img, itk::SmartPointer<TLabelImage> overlayImg);
  /*!
  	\fn Image::OverlayContour(itk::SmartPointer<TImage> img, itk::SmartPointer<TLabelImage> overlayImg)
  	\brief Generates an overlay of the images provided. Returns an RGB image of the result.
  */
  template<typename TLabelImage>
  static itk::SmartPointer< itk::Image<itk::RGBPixel<unsigned char>, TImage::ImageDimension> > OverlayContour(itk::SmartPointer<TImage> img, itk::SmartPointer<TLabelImage> overlayImg);
#endif // (ITK_REVIEW || ITK_VERSION_MAJOR > 3)

  //Thresholding
  /*!
  	\fn Image::ThresholdAboveImage(itk::SmartPointer<TImage> img, float outsideValue, float aboveValue)
  	\brief Generates an image with the intensities above a certain level thresholded (capped) to the outsideValue.
  */
  static itk::SmartPointer<TImage> ThresholdAboveImage(itk::SmartPointer<TImage> img, float outsideValue, float aboveValue);
  /*!
  	\fn Image::ThresholdBelowImage(itk::SmartPointer<TImage> img, float outsideValue, float aboveValue)
  	\brief Generates an image with the intensities below a certain level thresholded (capped) to the outsideValue.
  */
  static itk::SmartPointer<TImage> ThresholdBelowImage(itk::SmartPointer<TImage> img, float outsideValue, float aboveValue);
  /*!
  	\fn Image::ThresholdImage(itk::SmartPointer<TImage> img, float outsideValue, float belowValue, float aboveValue)
  	\brief Generates an image with the intensities below and above a certain level thresholded (capped) to the outsideValue.
  */
  static itk::SmartPointer<TImage> ThresholdImage(itk::SmartPointer<TImage> img, float outsideValue, float belowValue, float aboveValue);
  /*!
  	\fn Image::BinaryThresholdImage(itk::SmartPointer<TImage> img, float outsideValue, float insideValue, float belowValue, float aboveValue)
  	\brief Generates a binary image thresholded at the intensities below and above a certain level to the outsideValue.
  */
  template<typename TOutImage>
  static itk::SmartPointer<TOutImage> BinaryThresholdImage(itk::SmartPointer<TImage> img, float outsideValue, float insideValue, float belowValue, float aboveValue);
  /*!
  	\fn Image::OtsuThreshold(itk::SmartPointer<TImage> img, const int bins)
  	\brief Returns the Otsu threshold of an image of given the number of histogram bins.
  */
  static double OtsuThreshold(itk::SmartPointer<TImage> img, const int bins);
  /*!
  	\fn Image::OtsuThresholdImage(itk::SmartPointer<TImage> img, const int bins)
  	\brief Generates the Otsu threshold of an image of given the number of histogram bins.

    Uses binary threshold filter internally to threshold image above the Otsu threshold found.
  */
  template<typename TOutImage>
  static itk::SmartPointer<TOutImage> OtsuThresholdImage(itk::SmartPointer<TImage> img, const int bins);
  /*!
  	\fn Image::OtsuMultipleThresholdImage(itk::SmartPointer<TImage> img, const int bins, const int noOfLabels = 1)
  	\brief Generates the multiple Otsu threshold of an image of given the number of histogram bins.

    Uses binary threshold filter internally to threshold image above the multiple Otsu thresholds found.
  */
  template<typename TOutImage>
  static itk::SmartPointer<TOutImage> OtsuMultipleThresholdImage(itk::SmartPointer<TImage> img, const int bins, const int noOfLabels = 1);

  //Vector operations
  /*!
  	\fn Image::VectorMagnitude(itk::SmartPointer<TImage> img)
  	\brief Generates an image comprised of the magnitude of the vectors in the image.
  */
  template<typename TScalarImage>
  static itk::SmartPointer<TScalarImage> VectorMagnitude(itk::SmartPointer<TImage> img);
  //@}

  //Collection operations
  /**
  * \defgroup Image_Collection Image Collection Operations
  * \brief Members for operating of a collection of images (i.e. image batching)
  */
  //@{
  /*!
  	\fn Image::InformationCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Batch process images by output the image information for each image.
  */
  static void InformationCollection(std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
  	\fn Image::MatchHistogramCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TImage> refImage, const int bins = 128)
  	\brief Batch process images by matching the histograms of the images to the reference image provided.
  */
  static void MatchHistogramCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TImage> refImage, const int bins = 128);
  /*!
    \fn Image::CheckerboardCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TImage> refImage, const int squares = 10)
    \brief Batch process images by checkerboarding all of the images to the reference image provided.
  */
  static void CheckerboardCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TImage> refImage, const int squares = 10);
  /*!
  	\fn Image::RescaleIntensitiesCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float belowValue, float aboveValue)
  	\brief Batch process images by rescaling the intensities for each image.
  */
  static void RescaleIntensityCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float belowValue, float aboveValue);
  /*!
  	\fn Image::InvertIntensityCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Batch process images by inverting the intensities for each image.
  */
  static void InvertIntensityCollection(std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
  	\fn Image::FlipCollection(std::vector< typename itk::SmartPointer<TImage> > &images, bool xAxis = false, bool yAxis = true, bool zAxis = false, bool aboutOrigin = true)
  	\brief Batch process images by flipping each image along axis provided.
  */
  static void FlipCollection(std::vector< typename itk::SmartPointer<TImage> > &images, bool xAxis = false, bool yAxis = true, bool zAxis = false, bool aboutOrigin = true);
  /*!
  	\fn Image::AnisotropicDiffusionCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, const int iterations, float timestep = -1.0)
  	\brief Batch process images by smoothing each image using Gradient Anisotropic Diffusion.

  	Output has always has to be floating point.
  */
  template<typename TOutImage>
  static std::vector< typename itk::SmartPointer<TOutImage> > AnisotropicDiffusionCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, const int iterations, float timestep = -1.0);
  /*!
    \fn Image::BilateralCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float sigmaRange, float sigmaSpatial)
    \brief Batch process images by bilateral smoothing each image.
  */
  static void BilateralCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float sigmaRange, float sigmaSpatial);
  /*!
  	\fn Image::CastCollection(const std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Batch process images by casting each image to type provided in templates.
  */
  template<typename TOutImage>
  static std::vector< typename itk::SmartPointer<TOutImage> > CastCollection(const std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
  	\fn Image::MedianCollection(std::vector< typename itk::SmartPointer<TImage> > &images, const int radius)
  	\brief Batch process images by smoothing each image using the Median of the neighbourhood.
  */
  static void MedianCollection(std::vector< typename itk::SmartPointer<TImage> > &images, const int radius);
  /*!
  	\fn Image::GradientMagnitudeCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Batch process images by hightlighting the edges of each image using Gradient magnitude information of the pixels in the image.
  */
  static void GradientMagnitudeCollection(std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
  	\fn Image::LaplacianCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Batch process images by hightlighting the edges of each image using the Laplacian of the image.
  */
  static void LaplacianCollection(std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
  	\fn Image::DistanceMapCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, const bool binaryImage = false, const bool signedDistances = false, const bool computeInsideObject = false, const bool squaredDistance = false)
  	\brief Batch process images by computing the distance map of each image. Output should always be a float image.
  */
  template<typename TOutImage>
  static std::vector< typename itk::SmartPointer<TOutImage> > DistanceMapCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, const bool binaryImage = false, const bool signedDistances = false, const bool computeInsideObject = false, const bool squaredDistance = false);

  //Thresholding
  /*!
  	\fn Image::ThresholdAboveCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float aboveValue)
  	\brief Batch process images by thresholding each image from above.
  */
  static void ThresholdAboveCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float aboveValue);
  /*!
  	\fn Image::ThresholdBelowCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float belowValue)
  	\brief Batch process images by thresholding each image from below.
  */
  static void ThresholdBelowCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float belowValue);
  /*!
  	\fn Image::ThresholdCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float belowValue, float aboveValue)
  	\brief Batch process images by thresholding each image within band provided.
  */
  static void ThresholdCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float belowValue, float aboveValue);
  /*!
  	\fn Image::BinaryThresholdCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float insideValue, float belowValue, float aboveValue)
  	\brief Batch process images by binary thresholding each image within band and inside value provided.
  */
  template<typename TOutImage>
  static std::vector< typename itk::SmartPointer<TOutImage> > BinaryThresholdCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float insideValue, float belowValue, float aboveValue);
  /*!
  	\fn Image::OtsuMultipleThresholdCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, const int bins, const int noOfLabels)
  	\brief Batch process images by Otsu (histogram-based) thresholding each image having number of partitions/labels provided.
  */
  template<typename TOutImage>
  static std::vector< typename itk::SmartPointer<TOutImage> > OtsuMultipleThresholdCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, const int bins, const int noOfLabels);

  /*!
  \fn Image::SubsampleCollection(std::vector< typename itk::SmartPointer<TImage> > &images, unsigned factor = 2)
  \brief Batch process images by downsampling each image by factor provided.
  */
  static void SubsampleCollection(std::vector< typename itk::SmartPointer<TImage> > &images, unsigned factor = 2);
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
  /*!
  	\fn Image::MaskAndCropCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TMaskImage> maskImage, const size_t pixelPadding = 1)
  	\brief Batch process images by masking and cropping using the mask bounding box each image using the mask image provided.
  */
  template<class TMaskImage>
  static void MaskAndCropCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TMaskImage> maskImage, const size_t pixelPadding = 1);
#endif // (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
  /*!
  	\fn Image::MaskCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TMaskImage> maskImage)
  	\brief Batch process images by masking using the mask image provided.
  */
  template<class TMaskImage>
  static void MaskCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TMaskImage> maskImage);
  /*!
  	\fn Image::RelabelCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Batch process images by relabelling unconnected regions consecutatively.
  */
  static void RelabelCollection(std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
  	\fn Image::ResampleCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TRefImage> refImage)
  	\brief Batch process images by resampling images using the reference image provided.
  */
  template<class TRefImage>
  static void ResampleCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TRefImage> refImage);
  /*!
  	\fn Image::ResampleLabelCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TRefImage> refImage)
  	\brief Batch process labelled images by resampling images using the reference image provided. This uses the nearest neighbour interploator.
  */
  template<class TRefImage>
  static void ResampleLabelCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TRefImage> refImage);
  /*!
  	\fn Image::AddCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Batch process images by the adding the images together (pixel-wise).
  */
  static itk::SmartPointer<TImage> AddCollection(std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
  	\fn Image::SubtractCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Batch process images by the subtracting the images from the first image (pixel-wise).
  */
  static itk::SmartPointer<TImage> SubtractCollection(std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
  	\fn Image::DifferenceCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Same as SubtractCollection.
  */
  static inline itk::SmartPointer<TImage> DifferenceCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  {
    return SubtractCollection(images);
  }
  /*!
  	\fn Image::AverageCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
  	\brief Batch process images by the averaging the images (pixel-wise).
  */
  template<class TOutImage>
  static itk::SmartPointer<TOutImage> AverageCollection(std::vector< typename itk::SmartPointer<TImage> > &images);
  /*!
  	\fn Image::AverageVectorCollection(std::vector< typename itk::SmartPointer<TImage> > &images, int numberOfComponents)
  	\brief Batch process images by the averaging the vector images (pixel-wise).
  */
  static itk::SmartPointer<TImage> AverageVectorCollection(std::vector< typename itk::SmartPointer<TImage> > &images, int numberOfComponents);
  template<class TJoinedImage>
  itk::SmartPointer<TJoinedImage> JoinCollection(std::vector< typename itk::SmartPointer<TImage> > &images, const double origin, const double spacing);
#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
  /*!
  	\fn Image::MergeLabelledImages(std::vector< typename itk::SmartPointer<TImage> > &images, const size_t mergeType = 0)
  	\brief Batch process labelled images by the merging labels into one image.

    mergeType can be the following values:
    Type 0: Keep, Do its best to keep the label unchanged, but if a label is already used in a previous label map
    Type 1: Aggregate, If the same label is found several times in the label maps, the label objects with the same label are merged
    Type 2: Pack, Relabel all the label objects by order of processing
    Type 3: Strict, Keeps the labels unchanged and raises an exception if the same label is found in several images
  */
  static itk::SmartPointer<TImage> MergeLabelledImages(std::vector< typename itk::SmartPointer<TImage> > &images, const size_t mergeType = 0);
#endif // (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
  //@}

protected:
  itk::SmartPointer<TImage> CurrentImage; //!< Holds the current image in the pipeline
  itk::SmartPointer<TImage> PreviousImage; //!< Holds the previous image in the pipeline

private:

};

template<class TImage>
Image<TImage>::Image()
{

}

//template<class TImage>
//Image<TImage>::Image(itk::SmartPointer<TImage> image)
//{
//  SetInput(image);
//}

//template<class TImage>
//void Image<TImage>::SetInput(itk::SmartPointer<TImage> image)
//{
//  if(!CurrentImage)
//    CurrentImage = itk::SmartPointer<TImage>::New();
//  CurrentImage->DuplicateImage(image);
//}

//Arithmetic
template<class TImage>
vtkSmartPointer<vtkImageData> Image<TImage>::ConvertITKImageToVTKImage(itk::SmartPointer<TImage> img)
{
  typedef itk::ImageToVTKImageFilter<TImage> ConvertImageType;

  typename ConvertImageType::Pointer convertFilter = ConvertImageType::New();
  convertFilter->SetInput(img);
  convertFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    convertFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Converting ITK Image to VTK Image");
    PrintError(ex.GetDescription());
  }

  return convertFilter->GetOutput();
}

template<class TImage>
vtkSmartPointer<vtkImageData> Image<TImage>::ConvertITKVectorImageToVTKImage(itk::SmartPointer<TImage> img)
{
  typedef itk::Image<typename TImage::InternalPixelType, TImage::ImageDimension> TImageComponent;

  //Get the size etc. right of result image
  typename TImageComponent::Pointer componentImg = milx::Image<TImage>::ExtractComponent<TImageComponent>(img, 0);
  vtkSmartPointer<vtkImageData> fieldComponent = milx::Image<TImageComponent>::ConvertITKImageToVTKImage(componentImg);

  vtkSmartPointer<vtkImageData> field = vtkSmartPointer<vtkImageData>::New();
    field->SetDimensions(fieldComponent->GetDimensions());
    field->SetOrigin(fieldComponent->GetOrigin());
    field->SetSpacing(fieldComponent->GetSpacing());
  #if VTK_MAJOR_VERSION <= 5
    field->SetNumberOfScalarComponents(img->GetNumberOfComponentsPerPixel());
    field->SetScalarTypeToFloat();
    field->AllocateScalars();
  #else
    field->AllocateScalars(VTK_FLOAT, 3);
  #endif
  /*vtkSmartPointer<vtkImageAppendComponents> appendImage = vtkSmartPointer<vtkImageAppendComponents>::New(); //doesnt work!
    appendImage->DebugOn();
    appendImage->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates);
    appendImage->SetInput(fieldComponent);*/

  for(vtkIdType j = 0; j < img->GetNumberOfComponentsPerPixel(); j ++)
  {
    typename TImageComponent::Pointer componentImgTmp = milx::Image<TImage>::ExtractComponent<TImageComponent>(img, j);
    vtkSmartPointer<vtkImageData> fieldComponentTmp = milx::Image<TImageComponent>::ConvertITKImageToVTKImage(componentImgTmp);

    for(int k = 0; k < fieldComponent->GetDimensions()[0]; k ++)
      for(int l = 0; l < fieldComponent->GetDimensions()[1]; l ++)
        for(int m = 0; m < fieldComponent->GetDimensions()[2]; m ++)
          field->SetScalarComponentFromFloat(k, l, m, j, fieldComponentTmp->GetScalarComponentAsFloat(k, l, m, 0));
  }

  return field;
}

template<class TImage>
template<typename TPrecision>
itk::SmartPointer<TImage> Image<TImage>::ApplyOrientationToITKImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> refImage, const bool labelledImage, const bool flipY, const bool ignoreDirection)
{
  typename TImage::DirectionType direction = refImage->GetDirection();
  typename TImage::PointType origin = refImage->GetOrigin();

  typedef itk::InterpolateImageFunction<TImage, TPrecision> InterpolatorType;
  typename InterpolatorType::Pointer interpolator;
  if(labelledImage)
    interpolator = itk::NearestNeighborInterpolateImageFunction<TImage, TPrecision>::New();
  else
    interpolator = itk::LinearInterpolateImageFunction<TImage, TPrecision>::New();

  typename TImage::DirectionType identityMatrix;
  identityMatrix.SetIdentity();
//  typename TImage::PointType axesOrigin;
//  axesOrigin.Fill(0.0);

  //Remove origin and direction to get raw volume
  typedef itk::ChangeInformationImageFilter<TImage> ChangeFilterType;
  typename ChangeFilterType::Pointer stripInfo = ChangeFilterType::New();
    stripInfo->SetInput( img );
    stripInfo->ChangeDirectionOn();
    stripInfo->SetOutputDirection(identityMatrix);
//    stripInfo->ChangeOriginOn();
//    stripInfo->SetOutputOrigin(axesOrigin);

  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New(); //start with identity matrix
  matrix->Identity();
  for (int i = 0; i < 3; i ++)
    for (int k = 0; k < 3; k ++)
      matrix->SetElement(i,k, direction(i,k));
  matrix->Transpose(); //img volume to correct space, comment to go from space to img volume

  vtkSmartPointer<vtkMatrix4x4> transformMatrix = vtkSmartPointer<vtkMatrix4x4>::New(); //start with identity matrix
    transformMatrix->Identity();
    //Flip the image for VTK coordinate system
//    if(flipY)
//      transformMatrix->SetElement(1,1,-1); //flip

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Identity();
    transform->PostMultiply();
    transform->Concatenate(transformMatrix); //flip
    transform->Translate(-origin[0], -origin[1], -origin[2]); //remove origin
    transform->Concatenate(matrix); //direction
    transform->Translate(origin.GetDataPointer()); //add origin displacement
    transform->Update();
  vtkSmartPointer<vtkMatrix4x4> matrix2 = transform->GetMatrix();

  typename TImage::DirectionType orientMatrix;
  for (int i = 0; i < 3; i ++)
    for (int k = 0; k < 3; k ++)
      orientMatrix(i,k) = matrix2->GetElement(i,k);
  typename itk::Vector<TPrecision, milx::imgDimension> offset;
//  typename TImage::PointType new_origin;
  for (int i = 0; i < 3; i ++)
  {
//    new_origin[i] = matrix2->GetElement(i,3);
    offset[i] = matrix2->GetElement(i,3);
  }
  matrix2->Print(cout);

  typedef itk::ResampleImageFilter<TImage, TImage, TPrecision> ResampleImageFilterType;
  /*typename ResampleImageFilterType::Pointer stripInfo = ResampleImageFilterType::New();
    stripInfo->SetInput(img);
    stripInfo->SetDefaultPixelValue(0.0);
    stripInfo->SetSize( refImage->GetLargestPossibleRegion().GetSize() );
//    stripInfo->SetOutputOrigin(  axesOrigin );
    stripInfo->SetOutputSpacing( refImage->GetSpacing() );
//    stripInfo->SetOutputDirection( identityMatrix ); // 1 1 1 direction
    stripInfo->SetInterpolator(interpolator);
    stripInfo->AddObserver(itk::ProgressEvent(), ProgressUpdates);

  typedef itk::AffineTransform< TPrecision, milx::imgDimension > AffineTransformType;
  typedef typename AffineTransformType::Pointer AffineTransformPointer;
  AffineTransformPointer stripTransform = AffineTransformType::New();
  stripTransform->SetIdentity();
  typename TImage::DirectionType matrix = direction.GetInverse();
  stripTransform->SetMatrix(matrix);
  stripTransform->Translate(origin.GetVectorFromOrigin());
//  stripTransform->SetCenter(axesOrigin);
//  stripTransform->SetTranslation(origin.GetVectorFromOrigin());
  AffineTransformPointer stripInvTransform = AffineTransformType::New();
  stripInvTransform->GetInverse(stripTransform);

  stripInfo->SetTransform(stripInvTransform);*/

  typename TImage::SizeType inputSize = refImage->GetLargestPossibleRegion().GetSize();
  typename TImage::SizeType outputSize;
  typedef typename TImage::SizeType::SizeValueType SizeValueType;
  outputSize[0] = static_cast<SizeValueType>(inputSize[0] * 1);
  outputSize[1] = static_cast<SizeValueType>(inputSize[1] * 1);
  outputSize[2] = static_cast<SizeValueType>(inputSize[2] * 1);

  typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
    resample->SetInput(stripInfo->GetOutput());
//    resample->SetInput(img);
    resample->SetDefaultPixelValue(0.0);
    resample->SetSize( outputSize );
    resample->SetOutputStartIndex( refImage->GetLargestPossibleRegion().GetIndex() );
//    resample->SetOutputOrigin(  refImage->GetOrigin() );
    resample->SetOutputSpacing( refImage->GetSpacing() );
    resample->SetInterpolator(interpolator);
    resample->AddObserver(itk::ProgressEvent(), ProgressUpdates);

  ///Up cast transform to correct type
  typedef itk::AffineTransform< TPrecision, milx::imgDimension > AffineTransformType;
  typedef typename AffineTransformType::Pointer AffineTransformPointer;
  AffineTransformPointer affineTransform = AffineTransformType::New();
  affineTransform->SetIdentity();

  ///Transform
//  typename TImage::DirectionType matrix;
//  if(ignoreDirection)
//    matrix = identityMatrix;
//  else
//    matrix = direction.GetInverse();
//  if(flipY)
//    matrix(1,1) *= -1; //flip y
//  affineTransform->SetMatrix(matrix);
  affineTransform->SetMatrix(orientMatrix);

  affineTransform->SetOffset(offset);
//  affineTransform->SetCenter(origin);
//  affineTransform->SetTranslation(origin.GetVectorFromOrigin());

//  AffineTransformPointer invTransform = AffineTransformType::New();
//  invTransform->GetInverse(affineTransform);

  ///Origin
  AffineTransformPointer affineTransform2 = AffineTransformType::New();
  affineTransform2->SetIdentity();
//  typename TImage::DirectionType invOrientMatrix(orientMatrix.GetInverse());
//  affineTransform2->SetMatrix(invOrientMatrix);
  affineTransform2->SetMatrix(orientMatrix);
  //get zeroth point
  typename TImage::IndexType originIndex;
  originIndex.Fill(0);
  typename TImage::PointType old_origin;
  refImage->TransformIndexToPhysicalPoint(originIndex, old_origin);
  typename TImage::PointType new_origin = affineTransform2->TransformPoint(old_origin);
  resample->SetOutputOrigin(new_origin);
//  resample->SetOutputOrigin(origin);

  resample->SetTransform(affineTransform);

  try
  {
    resample->UpdateLargestPossibleRegion();
    resample->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed applying orientation to Image");
    PrintError(ex.GetDescription());
  }

  return resample->GetOutput();
//  stripInfo->Update();
//  return stripInfo->GetOutput();
}

#ifndef ITK_ONLY //Requires VTK
template<class TImage>
vtkSmartPointer<vtkImageData> Image<TImage>::ApplyOrientationToVTKImage(vtkSmartPointer<vtkImageData> img, itk::SmartPointer<TImage> refImage, vtkSmartPointer<vtkMatrix4x4> &transformMatrix, const bool labelledImage, const bool flipY)
{
  typename TImage::DirectionType direction = refImage->GetDirection();
  typename TImage::PointType origin = refImage->GetOrigin();

  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New(); //start with identity matrix
  matrix->Identity();
  for (int i = 0; i < 3; i ++)
    for (int k = 0; k < 3; k ++)
      matrix->SetElement(i,k, direction(i,k));
  matrix->Transpose(); //img volume to correct space, comment to go from space to img volume

  if(!transformMatrix)
    transformMatrix = vtkSmartPointer<vtkMatrix4x4>::New(); //start with identity matrix
  transformMatrix->Identity();
  //Flip the image for VTK coordinate system
  if(flipY)
    transformMatrix->SetElement(1,1,-1); //flip

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->Identity();
  transform->PostMultiply();
  transform->Concatenate(transformMatrix); //flip
  transform->Translate(-origin[0], -origin[1], -origin[2]); //remove origin
  transform->Concatenate(matrix); //direction
  transform->Translate(origin.GetDataPointer()); //add origin displacement

  vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
#if VTK_MAJOR_VERSION <= 5
  reslice->SetInput(img);
#else
  reslice->SetInputData(img);
#endif
  reslice->AddObserver(vtkCommand::ProgressEvent, VTKProgressUpdates);
  reslice->AutoCropOutputOn();
  reslice->TransformInputSamplingOn();
  reslice->SetOutputDimensionality(milx::imgDimension);
  reslice->SetResliceAxes(transform->GetMatrix());

  if(!labelledImage)
      reslice->SetInterpolationModeToLinear(); //reduce artefacts for normal images
  else
  {
      reslice->SetInterpolationModeToNearestNeighbor(); //reduce artefacts for labels
      PrintDebug("Using NN Interpolator for Image Orientation.");
//  #if VTK_MAJOR_VERSION > 5
//      reslice->SetOutputScalarType(VTK_CHAR);
//  #endif
  }
  reslice->Update();

  std::string resliceType = reslice->GetOutput()->GetScalarTypeAsString();
  PrintDebug("Orientated Image Type as " + resliceType);

  return reslice->GetOutput();
}
#endif

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::ConvertVTKImageToITKImage(vtkSmartPointer<vtkImageData> img)
{
  typedef itk::VTKImageToImageFilter<TImage> ConvertImageType;

  typename ConvertImageType::Pointer convertFilter = ConvertImageType::New();
  convertFilter->SetInput(img);
  convertFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    convertFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Converting VTK Image to ITK Image");
    PrintError(ex.GetDescription());
  }

  itk::SmartPointer<TImage> tmpImage = DuplicateImage(convertFilter->GetOutput());

  return tmpImage;
}

template<class TImage>
template<typename TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::CastImage(itk::SmartPointer<TImage> img)
{
  typedef itk::CastImageFilter<TImage, TOutImage> CastImageType;

  typename CastImageType::Pointer castFilter = CastImageType::New();
  castFilter->SetInput(img);
  castFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    castFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Casting");
    PrintError(ex.GetDescription());
  }

  return castFilter->GetOutput();
}

/*template<class TImage>
template<typename TOutImage>
void Image<TImage>::DeepCopy(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToCopyFrom, TImage::RegionType region)
{
  itk::ImageRegionConstIterator<TImage> inputIterator(imgToCopyFrom, region);
  itk::ImageRegionIterator<TImage> outputIterator(img, region);
  while(!inputIterator.IsAtEnd())
    {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
    }
}*/

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::DuplicateImage(itk::SmartPointer<TImage> img)
{
  typedef itk::ImageDuplicator<TImage> DuplicateType;

  typename DuplicateType::Pointer duplicator = DuplicateType::New();
  duplicator->SetInputImage(img);
  duplicator->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    duplicator->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Generating Duplicate Image");
    PrintError(ex.GetDescription());
  }

  return duplicator->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::DuplicateImage(const TImage *img)
{
  typedef itk::ImageDuplicator<TImage> DuplicateType;

  typename DuplicateType::Pointer duplicator = DuplicateType::New();
  duplicator->SetInputImage(img);
  duplicator->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    duplicator->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Generating Duplicate Image");
    PrintError(ex.GetDescription());
  }

  return duplicator->GetOutput();
}

template<class TImage>
template<typename TVector>
itk::SmartPointer<TImage> Image<TImage>::ImportVectorToImage(vnl_vector<TVector> &vec, typename TImage::SizeType size, itk::SmartPointer<TImage> image)
{
  //Pass by reference is necessary here.
  typedef itk::ImportImageFilter<typename TImage::PixelType, TImage::ImageDimension> ImportFilterType;
  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();

  typename ImportFilterType::IndexType start;
  start.Fill( 0 );

  typename ImportFilterType::RegionType region;
  region.SetIndex(start);
  region.SetSize(size);

  importFilter->SetRegion( region );
  if(image)
  {
    importFilter->SetOrigin( image->GetOrigin() );
    importFilter->SetSpacing( image->GetSpacing() );
    importFilter->SetDirection( image->GetDirection() );
  }

  const bool importImageFilterWillOwnTheBuffer = false;
  importFilter->SetImportPointer( vec.data_block(), vec.size(), importImageFilterWillOwnTheBuffer );
  importFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    importFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Converting vector to ITK Image");
    PrintError(ex.GetDescription());
  }

  return importFilter->GetOutput();
}

template<class TImage>
template<typename TMatrix>
itk::SmartPointer<TImage> Image<TImage>::ImportMatrixToImage(vnl_matrix<TMatrix> &matrix, itk::SmartPointer<TImage> image)
{
  //Pass by reference is necessary here.
  typedef itk::ImportImageFilter<typename TImage::PixelType, TImage::ImageDimension> ImportFilterType;
  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();

  typename ImportFilterType::SizeType  size;
  size[0]  = matrix.rows();  // size along X
  size[1]  = matrix.cols();  // size along Y
  size[2]  = 1;  // size along Z

  typename ImportFilterType::IndexType start;
  start.Fill( 0 );

  typename ImportFilterType::RegionType region;
  region.SetIndex(start);
  region.SetSize(size);

  importFilter->SetRegion( region );
  if(image)
  {
    importFilter->SetOrigin( image->GetOrigin() );
    importFilter->SetSpacing( image->GetSpacing() );
    importFilter->SetDirection( image->GetDirection() );
  }

  const bool importImageFilterWillOwnTheBuffer = false;
  importFilter->SetImportPointer( matrix.data_block(), matrix.size(), importImageFilterWillOwnTheBuffer );
  importFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    importFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Converting matrix to ITK Image");
    PrintError(ex.GetDescription());
  }

  return importFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::BlankImage(typename TImage::PixelType value, typename TImage::SizeType imgSize)
{
  typename TImage::IndexType blankStart;
  blankStart.Fill(0);

  typename TImage::RegionType blankRegion;
  blankRegion.SetIndex(blankStart);
  blankRegion.SetSize(imgSize);

  //create image
  typename TImage::Pointer img = TImage::New();
  img->SetRegions(blankRegion);
  img->Allocate();
  img->FillBuffer(value);

  return img;
}

//Geometry
template<class TImage>
template<typename TSubImage>
itk::SmartPointer<TSubImage> Image<TImage>::ExtractSubImage(itk::SmartPointer<TImage> img, typename TImage::RegionType imgRegion)
{
  typedef itk::ExtractImageFilter<TImage, TSubImage> FilterType;
  typename FilterType::Pointer extractor = FilterType::New();
  extractor->SetExtractionRegion(imgRegion);
  extractor->SetInput(img);
#if ITK_VERSION_MAJOR >= 4
  extractor->SetDirectionCollapseToIdentity(); // This is required.
#endif
  extractor->AddObserver(itk::ProgressEvent(), ProgressUpdates);

  try
  {
    extractor->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Extracting Sub-Image");
    PrintError(ex.GetDescription());
  }

  return extractor->GetOutput();
}

template<class TImage>
template<typename TImageSlice>
itk::SmartPointer<TImageSlice> Image<TImage>::ExtractSlice(itk::SmartPointer<TImage> img, int *extent)
{
  typename TImage::IndexType desiredStart; ///Pick the hyper plane
    desiredStart[0] = extent[0];
    desiredStart[1] = extent[2];
    desiredStart[2] = extent[4];
  milx::PrintDebug("Slice Start Indices: " + milx::NumberToString(desiredStart[0]) + ", " + milx::NumberToString(desiredStart[1]) + ", " + milx::NumberToString(desiredStart[2]));

  typename TImage::SizeType desiredSize;
    desiredSize[0] = extent[1]-extent[0];
    desiredSize[1] = extent[3]-extent[2];
    desiredSize[2] = extent[5]-extent[4];

  const size_t sliceDim = TImageSlice::ImageDimension;
  const size_t imgDim = TImage::ImageDimension;

  if(sliceDim == imgDim) //then just slicing data and not compacting a diemnsion
  {
    for(size_t j = 0; j < imgDim; j ++)
    {
      if(desiredSize[j] == 0) //dont compact
        desiredSize[j] += 1;
    }
  }

  milx::PrintDebug("Slice Size: " + milx::NumberToString(desiredSize[0]) + ", " + milx::NumberToString(desiredSize[1]) + ", " + milx::NumberToString(desiredSize[2]));

  typename TImage::RegionType desiredRegion(desiredStart, desiredSize);

  return ExtractSubImage<TImageSlice>(img, desiredRegion);
}

template<class TImage>
template<typename TImageComponent>
itk::SmartPointer<TImageComponent> Image<TImage>::ExtractComponent(itk::SmartPointer<TImage> img, int component)
{
  typedef itk::VectorIndexSelectionCastImageFilter<TImage, TImageComponent> IndexSelectionType;
  typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(component);
    indexSelectionFilter->SetInput(img);

  try
  {
    indexSelectionFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Extracting Component");
    PrintError(ex.GetDescription());
  }

  return indexSelectionFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::ResizeImage(itk::SmartPointer<TImage> img, typename TImage::SizeType imgSize)
{
  typedef itk::IdentityTransform<double, imgDimension> TransformType;
  typedef itk::ResampleImageFilter<TImage, TImage> ResampleImageFilterType;
  typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
  resample->SetInput(img);
  resample->SetSize(imgSize);
  resample->SetOutputSpacing(img->GetSpacing());
  resample->SetTransform(TransformType::New());
  resample->UpdateLargestPossibleRegion();
  resample->AddObserver(itk::ProgressEvent(), ProgressUpdates);

  try
  {
    resample->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Resizing");
    PrintError(ex.GetDescription());
  }

  return resample->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::ResizeImage(itk::SmartPointer<TImage> img, typename TImage::SizeType imgSize, typename TImage::SpacingType outputSpacing)
{
  typedef itk::IdentityTransform<double, imgDimension> TransformType;
  typedef itk::ResampleImageFilter<TImage, TImage> ResampleImageFilterType;
  typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
  resample->SetInput(img);
  resample->SetSize(imgSize);
  resample->SetOutputSpacing(outputSpacing);
  resample->SetTransform(TransformType::New());
  resample->UpdateLargestPossibleRegion();
  resample->AddObserver(itk::ProgressEvent(), ProgressUpdates);

  try
  {
    resample->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Resizing");
    PrintError(ex.GetDescription());
  }

  return resample->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::ResizeImage(itk::SmartPointer<TImage> img, typename TImage::SizeType imgSize, typename TImage::SpacingType outputSpacing, typename TImage::PointType outputOrigin, typename TImage::DirectionType outputDirection)
{
  typedef itk::IdentityTransform<double, imgDimension> TransformType;
  typedef itk::ResampleImageFilter<TImage, TImage> ResampleImageFilterType;
  typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
  resample->SetInput(img);
  resample->SetSize(imgSize);
  resample->SetOutputSpacing(outputSpacing);
  resample->SetOutputOrigin(outputOrigin);
  resample->SetOutputDirection(outputDirection);
  resample->SetTransform(TransformType::New());
  resample->UpdateLargestPossibleRegion();
  resample->AddObserver(itk::ProgressEvent(), ProgressUpdates);

  try
  {
    resample->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Resizing");
    PrintError(ex.GetDescription());
  }

  return resample->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::SubsampleImage(itk::SmartPointer<TImage> img, typename TImage::SizeType factors)
{
  typedef itk::ShrinkImageFilter<TImage, TImage> SubSampleImageFilterType;
  typename SubSampleImageFilterType::Pointer subsample = SubSampleImageFilterType::New();
  subsample->SetInput(img);
  for(size_t j = 0; j < TImage::ImageDimension; j ++)
    subsample->SetShrinkFactor(j, factors[j]); // shrink the first dimension by a factor of 2 (i.e. 100 gets changed to 50)
  subsample->AddObserver(itk::ProgressEvent(), ProgressUpdates);

  try
  {
    subsample->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Subsampling");
    PrintError(ex.GetDescription());
  }

  return subsample->GetOutput();
}

//Transform
template<class TImage>
template<typename TOutImage, typename TTransform, typename TPrecision>
itk::SmartPointer<TOutImage> Image<TImage>::TransformImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TOutImage> refImg, itk::SmartPointer<TTransform> transf, const bool inverse, const int interp)
{
  typedef itk::InterpolateImageFunction<TImage, TPrecision> InterpolatorType;
  typename InterpolatorType::Pointer interpolator;
  if(interp == 0)
  {
    interpolator = itk::NearestNeighborInterpolateImageFunction<TImage, TPrecision>::New();
    PrintDebug("Using nearest neighbour interpolation for resampling.");
  }
  else if(interp == 1)
  {
    interpolator = itk::LinearInterpolateImageFunction<TImage, TPrecision>::New();
    PrintDebug("Using linear interpolation for resampling.");
  }
  else
  {
    typedef itk::BSplineInterpolateImageFunction<TImage, TPrecision, TPrecision> BSplineInterpolatorType;
    typename BSplineInterpolatorType::Pointer bsplineInterpolator = BSplineInterpolatorType::New();
    bsplineInterpolator->SetSplineOrder(3);
    interpolator = bsplineInterpolator;
    PrintDebug("Using BSpline interpolation for resampling.");
  }

  typedef itk::ResampleImageFilter<TImage, TOutImage, TPrecision> ResampleImageFilterType;
  typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
  resample->SetInput(img);
  resample->SetDefaultPixelValue(0.0);
  resample->UseReferenceImageOn();
  resample->SetReferenceImage(refImg);
//  resample->SetSize( refImg->GetLargestPossibleRegion().GetSize() );
//  resample->SetOutputOrigin(  refImg->GetOrigin() );
//  resample->SetOutputSpacing( refImg->GetSpacing() );
//  resample->SetOutputDirection( refImg->GetDirection() );
  resample->SetInterpolator(interpolator);
  resample->AddObserver(itk::ProgressEvent(), ProgressUpdates);

  ///Up cast transform to correct type
  typedef itk::AffineTransform< TPrecision, milx::imgDimension > AffineTransformType;
  typedef typename AffineTransformType::Pointer AffineTransformPointer;
  typedef itk::VersorRigid3DTransform<TPrecision> VersorRigidTransformType;
  typedef typename VersorRigidTransformType::Pointer VersorRigidTransformPointer;
#if ITK_VERSION_MAJOR > 3 //New() member issues, see ITK4 tests for detail
  typedef itkv3::Rigid3DTransform<TPrecision> RigidTransformType;
#else
  typedef itk::Rigid3DTransform<TPrecision> RigidTransformType;
#endif
  typedef typename RigidTransformType::Pointer RigidTransformPointer;
  typedef itk::CenteredEuler3DTransform<TPrecision> CenteredEuler3DTransformType;
  typedef typename CenteredEuler3DTransformType::Pointer CenteredEuler3DTransformPointer;
  typedef itk::Euler3DTransform<TPrecision> Euler3DTransformType;
  typedef typename Euler3DTransformType::Pointer Euler3DTransformPointer;

  AffineTransformPointer affineTransform;
  VersorRigidTransformPointer versorRigidTransform;
  RigidTransformPointer rigidTransform;
  CenteredEuler3DTransformPointer centeredEulerTransform;
  Euler3DTransformPointer euler3DTransform;

  if( !strcmp(transf->GetNameOfClass(),"AffineTransform") )
  {
    AffineTransformPointer affine_read = static_cast<AffineTransformType*>(transf.GetPointer());
    affineTransform = dynamic_cast< AffineTransformType * >( affine_read.GetPointer() );
    if( affineTransform )
    {
      PrintInfo("Successfully Read Affine Transform file.");
      AffineTransformPointer affineInvTransform;
      if(inverse)
      {
        affineInvTransform = AffineTransformType::New();
        affineTransform->GetInverse(affineInvTransform);
      }
      else
        affineInvTransform = affineTransform;
      resample->SetTransform(affineInvTransform);
    }
  }
  else if( !strcmp(transf->GetNameOfClass(),"VersorRigid3DTransform") )
  {
    VersorRigidTransformPointer rigid_read = static_cast<VersorRigidTransformType*>(transf.GetPointer());
    versorRigidTransform = dynamic_cast< VersorRigidTransformType * >( rigid_read.GetPointer() );
    if( versorRigidTransform )
    {
      PrintInfo("Successfully Read Versor Rigid Transform file.");
      VersorRigidTransformPointer versorRigidInvTransform;
      if(inverse)
      {
        versorRigidInvTransform = VersorRigidTransformType::New();
        versorRigidTransform->GetInverse(versorRigidInvTransform);
      }
      else
        versorRigidInvTransform = versorRigidTransform;
      resample->SetTransform(versorRigidInvTransform);
    }
  }
  else if( !strcmp(transf->GetNameOfClass(),"Rigid3DTransform") )
  {
    RigidTransformPointer rigid_read = static_cast<RigidTransformType*>(transf.GetPointer());
    rigidTransform = dynamic_cast< RigidTransformType * >( rigid_read.GetPointer() );
    if( rigidTransform )
    {
      PrintInfo("Successfully Read Rigid Transform file.");
      RigidTransformPointer rigidInvTransform;
      if(inverse)
      {
        rigidInvTransform = RigidTransformType::New();
        rigidTransform->GetInverse(rigidInvTransform);
      }
      else
        rigidInvTransform = rigidTransform;
      resample->SetTransform(rigidInvTransform);
    }
  }
  else if( !strcmp(transf->GetNameOfClass(),"CenteredEuler3DTransform") )
  {
    CenteredEuler3DTransformPointer rigid_read = static_cast<CenteredEuler3DTransformType*>(transf.GetPointer());
    centeredEulerTransform = dynamic_cast< CenteredEuler3DTransformType * >( rigid_read.GetPointer() );
    if( centeredEulerTransform )
    {
      PrintInfo("Successfully Read Centered Euler 3D Transform file.");
      CenteredEuler3DTransformPointer centeredEulerInvTransform;
      if(inverse)
      {
        centeredEulerInvTransform = CenteredEuler3DTransformType::New();
        centeredEulerTransform->GetInverse(centeredEulerInvTransform);
      }
      else
        centeredEulerInvTransform = centeredEulerTransform;
      resample->SetTransform(centeredEulerInvTransform);
    }
  }
  else if( !strcmp(transf->GetNameOfClass(),"Euler3DTransform") )
  {
    Euler3DTransformPointer rigid_read = static_cast<Euler3DTransformType*>(transf.GetPointer());
    euler3DTransform = dynamic_cast< Euler3DTransformType * >( rigid_read.GetPointer() );
    if( euler3DTransform )
    {
      PrintInfo("Successfully Read Euler 3D Transform file.");
      Euler3DTransformPointer euler3DInvTransform;
      if(inverse)
      {
        euler3DInvTransform = Euler3DTransformType::New();
        euler3DTransform->GetInverse(euler3DInvTransform);
      }
      else
        euler3DInvTransform = euler3DTransform;
      resample->SetTransform(euler3DInvTransform);
    }
  }
  else
  {
      PrintError("milxImage: Transform Input is not of known type");
  }

  try
  {
    resample->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Transforming Image");
    PrintError(ex.GetDescription());
  }

  return resample->GetOutput();
}

template<class TImage>
template<typename TOutImage, typename TTransform, typename TPrecision>
itk::SmartPointer<TOutImage> Image<TImage>::TransformImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TTransform> transf, const bool inverse, const int interp)
{
  typedef itk::InterpolateImageFunction<TImage, TPrecision> InterpolatorType;
  typename InterpolatorType::Pointer interpolator;
  if(interp == 0)
  {
    interpolator = itk::NearestNeighborInterpolateImageFunction<TImage, TPrecision>::New();
    PrintDebug("Using nearest neighbour interpolation for resampling.");
  }
  else if(interp == 1)
  {
    interpolator = itk::LinearInterpolateImageFunction<TImage, TPrecision>::New();
    PrintDebug("Using linear interpolation for resampling.");
  }
  else
  {
    typedef itk::BSplineInterpolateImageFunction<TImage, TPrecision, TPrecision> BSplineInterpolatorType;
    typename BSplineInterpolatorType::Pointer bsplineInterpolator = BSplineInterpolatorType::New();
    bsplineInterpolator->SetSplineOrder(3);
    interpolator = bsplineInterpolator;
    PrintDebug("Using BSpline interpolation for resampling.");
  }

  //transform origin
  typename TImage::PointType transformedOrigin = transf->TransformPoint(img->GetOrigin());

  typedef itk::ResampleImageFilter<TImage, TOutImage, TPrecision> ResampleImageFilterType;
  typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();

  ///Up cast transform to correct type
  typedef itk::AffineTransform< TPrecision, milx::imgDimension > AffineTransformType;
  typedef typename AffineTransformType::Pointer AffineTransformPointer;
  typedef itk::VersorRigid3DTransform<TPrecision> VersorRigidTransformType;
  typedef typename VersorRigidTransformType::Pointer VersorRigidTransformPointer;
#if ITK_VERSION_MAJOR > 3 //New() member issues, see ITK4 tests for detail
  typedef itkv3::Rigid3DTransform<TPrecision> RigidTransformType;
#else
  typedef itk::Rigid3DTransform<TPrecision> RigidTransformType;
#endif
  typedef typename RigidTransformType::Pointer RigidTransformPointer;
  typedef itk::CenteredEuler3DTransform<TPrecision> CenteredEuler3DTransformType;
  typedef typename CenteredEuler3DTransformType::Pointer CenteredEuler3DTransformPointer;
  typedef itk::Euler3DTransform<TPrecision> Euler3DTransformType;
  typedef typename Euler3DTransformType::Pointer Euler3DTransformPointer;

  AffineTransformPointer affineTransform;
  VersorRigidTransformPointer versorRigidTransform;
  RigidTransformPointer rigidTransform;
  CenteredEuler3DTransformPointer centeredEulerTransform;
  Euler3DTransformPointer euler3DTransform;

  if( !strcmp(transf->GetNameOfClass(),"AffineTransform") )
  {
    AffineTransformPointer affine_read = static_cast<AffineTransformType*>(transf.GetPointer());
    affineTransform = dynamic_cast< AffineTransformType * >( affine_read.GetPointer() );
    if( affineTransform )
    {
      PrintInfo("Successfully Read Affine Transform file.");
      AffineTransformPointer affineInvTransform;
      if(inverse)
      {
        affineInvTransform = AffineTransformType::New();
        affineTransform->GetInverse(affineInvTransform);
      }
      else
        affineInvTransform = affineTransform;
      resample->SetTransform(affineInvTransform);
      transformedOrigin = affineInvTransform->TransformPoint(img->GetOrigin());
    }
  }
  else if( !strcmp(transf->GetNameOfClass(),"VersorRigid3DTransform") )
  {
    VersorRigidTransformPointer rigid_read = static_cast<VersorRigidTransformType*>(transf.GetPointer());
    versorRigidTransform = dynamic_cast< VersorRigidTransformType * >( rigid_read.GetPointer() );
    if( versorRigidTransform )
    {
      PrintInfo("Successfully Read Versor Rigid Transform file.");
      VersorRigidTransformPointer versorRigidInvTransform;
      if(inverse)
      {
        versorRigidInvTransform = VersorRigidTransformType::New();
        versorRigidTransform->GetInverse(versorRigidInvTransform);
      }
      else
        versorRigidInvTransform = versorRigidTransform;
      resample->SetTransform(versorRigidInvTransform);
      transformedOrigin = versorRigidInvTransform->TransformPoint(img->GetOrigin());
    }
  }
  else if( !strcmp(transf->GetNameOfClass(),"Rigid3DTransform") )
  {
    RigidTransformPointer rigid_read = static_cast<RigidTransformType*>(transf.GetPointer());
    rigidTransform = dynamic_cast< RigidTransformType * >( rigid_read.GetPointer() );
    if( rigidTransform )
    {
      PrintInfo("Successfully Read Rigid Transform file.");
      RigidTransformPointer rigidInvTransform;
      if(inverse)
      {
        rigidInvTransform = RigidTransformType::New();
        rigidTransform->GetInverse(rigidInvTransform);
      }
      else
        rigidInvTransform = rigidTransform;
      resample->SetTransform(rigidInvTransform);
      transformedOrigin = rigidInvTransform->TransformPoint(img->GetOrigin());
    }
  }
  else if( !strcmp(transf->GetNameOfClass(),"CenteredEuler3DTransform") )
  {
    CenteredEuler3DTransformPointer rigid_read = static_cast<CenteredEuler3DTransformType*>(transf.GetPointer());
    centeredEulerTransform = dynamic_cast< CenteredEuler3DTransformType * >( rigid_read.GetPointer() );
    if( centeredEulerTransform )
    {
      PrintInfo("Successfully Read Centered Euler 3D Transform file.");
      CenteredEuler3DTransformPointer centeredEulerInvTransform;
      if(inverse)
      {
        centeredEulerInvTransform = CenteredEuler3DTransformType::New();
        centeredEulerTransform->GetInverse(centeredEulerInvTransform);
      }
      else
        centeredEulerInvTransform = centeredEulerTransform;
      resample->SetTransform(centeredEulerInvTransform);
      transformedOrigin = centeredEulerInvTransform->TransformPoint(img->GetOrigin());
    }
  }
  else if( !strcmp(transf->GetNameOfClass(),"Euler3DTransform") )
  {
    Euler3DTransformPointer rigid_read = static_cast<Euler3DTransformType*>(transf.GetPointer());
    euler3DTransform = dynamic_cast< Euler3DTransformType * >( rigid_read.GetPointer() );
    if( euler3DTransform )
    {
      PrintInfo("Successfully Read Euler 3D Transform file.");
      Euler3DTransformPointer euler3DInvTransform;
      if(inverse)
      {
        euler3DInvTransform = Euler3DTransformType::New();
        euler3DTransform->GetInverse(euler3DInvTransform);
      }
      else
        euler3DInvTransform = euler3DTransform;
      resample->SetTransform(euler3DInvTransform);
      transformedOrigin = euler3DInvTransform->TransformPoint(img->GetOrigin());
    }
  }
  else
  {
      PrintError("milxImage: Transform Input is not of known type");
  }

  resample->SetInput(img);
  resample->SetDefaultPixelValue(0.0);
  resample->SetSize( img->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin( transformedOrigin );
  resample->SetOutputSpacing( img->GetSpacing() );
  resample->SetOutputDirection( img->GetDirection() );
  resample->SetInterpolator(interpolator);
  resample->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    resample->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Transforming Image");
    PrintError(ex.GetDescription());
  }

  return resample->GetOutput();
}

template<class TImage>
template<typename TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::ResampleImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TOutImage> refImg, const bool linearInterp)
{
//  typedef itk::IdentityTransform<float, imgDimension> TransformType;
  typedef itk::InterpolateImageFunction<TImage, float> InterpolatorType;
  typename InterpolatorType::Pointer interpolator;
  if(linearInterp)
  {
    interpolator = itk::LinearInterpolateImageFunction<TImage, float>::New();
    PrintDebug("Using linear interpolation for resampling.");
  }
  else
  {
    typedef itk::BSplineInterpolateImageFunction<TImage, float, float> BSplineInterpolatorType;
    typename BSplineInterpolatorType::Pointer bsplineInterpolator = BSplineInterpolatorType::New();
    bsplineInterpolator->SetSplineOrder(3);
    interpolator = bsplineInterpolator;
    PrintDebug("Using BSpline interpolation for resampling.");
  }

  typedef itk::ResampleImageFilter<TImage, TOutImage, float> ResampleImageFilterType;
  typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
  resample->SetInput(img);
  resample->SetDefaultPixelValue(0.0);
  resample->UseReferenceImageOn();
  resample->SetReferenceImage(refImg);
//  resample->SetOutputSpacing(refImg->GetSpacing());
//  resample->SetOutputOrigin(refImg->GetOrigin());
//  resample->SetOutputDirection(refImg->GetDirection());
//  resample->SetTransform(TransformType::New());
  resample->SetInterpolator(interpolator);
  resample->AddObserver(itk::ProgressEvent(), ProgressUpdates);
//  if(matchSize)
//  {
//    PrintDebug("Resampling with size changes.");
//    resample->SetSize(refImg->GetLargestPossibleRegion().GetSize());
//    resample->SetOutputStartIndex(refImg->GetLargestPossibleRegion().GetIndex());
//  }

  try
  {
    resample->UpdateLargestPossibleRegion();
    resample->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Resampling");
    PrintError(ex.GetDescription());
  }

  return resample->GetOutput();
}

template<class TImage>
template<typename TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::ResampleLabel(itk::SmartPointer<TImage> img, itk::SmartPointer<TOutImage> refImg)
{
//  typedef itk::IdentityTransform<float, imgDimension> TransformType;
  typedef itk::NearestNeighborInterpolateImageFunction<TImage, float> NNInterpolatorType;
  typename NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();
  PrintDebug("Using Nearest Neighbour interpolation for resampling.");

  typedef itk::ResampleImageFilter<TImage, TOutImage, float> ResampleImageFilterType;
  typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
  resample->SetInput(img);
  resample->SetDefaultPixelValue(0.0);
  resample->UseReferenceImageOn();
  resample->SetReferenceImage(refImg);
//  resample->SetOutputSpacing(refImg->GetSpacing());
//  resample->SetOutputOrigin(refImg->GetOrigin());
//  resample->SetOutputDirection(refImg->GetDirection());
//  resample->SetTransform(TransformType::New());
  resample->SetInterpolator(interpolator);
  resample->AddObserver(itk::ProgressEvent(), ProgressUpdates);
//  if(matchSize)
//  {
//    PrintDebug("Resampling with size changes.");
//    resample->SetSize(refImg->GetLargestPossibleRegion().GetSize());
//    resample->SetOutputStartIndex(refImg->GetLargestPossibleRegion().GetIndex());
//  }

  try
  {
    resample->UpdateLargestPossibleRegion();
    resample->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Resampling");
    PrintError(ex.GetDescription());
  }

  return resample->GetOutput();
}

//Arithmetic
template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::AddImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2)
{
  typedef itk::AddImageFilter<TImage, TImage> AddImageType;

  typename AddImageType::Pointer addFilter = AddImageType::New();
  addFilter->SetInput1(img1);
  addFilter->SetInput2(img2);
  addFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    addFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Addition");
    PrintError(ex.GetDescription());
  }

  return addFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::SubtractImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2)
{
  typedef itk::SubtractImageFilter<TImage, TImage> SubtractImageType;

  typename SubtractImageType::Pointer subtractFilter = SubtractImageType::New();
  subtractFilter->SetInput1(img1);
  subtractFilter->SetInput2(img2);
  subtractFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    subtractFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Subtraction");
    PrintError(ex.GetDescription());
  }

  return subtractFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::MultiplyImages(itk::SmartPointer<TImage> img1, itk::SmartPointer<TImage> img2)
{
  typedef itk::MultiplyImageFilter<TImage, TImage> MultiplyImageType;

  typename MultiplyImageType::Pointer multiplyFilter = MultiplyImageType::New();
  multiplyFilter->SetInput1(img1);
  multiplyFilter->SetInput2(img2);
  multiplyFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    multiplyFilter->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    PrintError("Failed Multiplication");
    PrintError(ex.GetDescription());
  }

  return multiplyFilter->GetOutput();
}

template<class TImage>
template<class TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::ScaleImage(itk::SmartPointer<TImage> img, float scaling)
{
  typedef itk::ShiftScaleImageFilter<TImage, TOutImage> ScaleImageType;

  typename ScaleImageType::Pointer scaleFilter = ScaleImageType::New();
  scaleFilter->SetInput(img);
  scaleFilter->SetScale(scaling);
  scaleFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    scaleFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Scaling");
    PrintError(ex.GetDescription());
  }

  return scaleFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::ScaleVectorImage(itk::SmartPointer<TImage> img, float scaling, int numberOfComponents)
{
  typedef itk::VectorShiftScaleImageFilter<TImage, TImage> ScaleImageType;

  typename ScaleImageType::Pointer scaleFilter = ScaleImageType::New();
  scaleFilter->SetInput(img);
  scaleFilter->SetScale(scaling);
  scaleFilter->SetDimension(numberOfComponents);
  scaleFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    scaleFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Scaling");
    PrintError(ex.GetDescription());
  }

  return scaleFilter->GetOutput();
}

#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::ConvolveImages(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> kernelImg)
{
#if ITK_VERSION_MAJOR > 3
  typedef itk::FFTConvolutionImageFilter<TImage, TImage> ConvolveImageType;
#else
  typedef itk::ConvolutionImageFilter<TImage, TImage> ConvolveImageType;
#endif

  typename ConvolveImageType::Pointer convolveFilter = ConvolveImageType::New();
  convolveFilter->SetInput(img);
#if ITK_VERSION_MAJOR > 3
  convolveFilter->SetKernelImage(kernelImg);
#else
  convolveFilter->SetImageKernelInput(kernelImg);
#endif
  convolveFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    convolveFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Convolution");
    PrintError(ex.GetDescription());
  }

  return convolveFilter->GetOutput();
}
#endif //(ITK_REVIEW || ITK_VERSION_MAJOR > 3)

//Info
template<class TImage>
void Image<TImage>::Information(itk::SmartPointer<TImage> img)
{
  typename TImage::DirectionType direction = img->GetDirection();
  typename TImage::PointType origin = img->GetOrigin();
  typename TImage::SpacingType spacing = img->GetSpacing();
  typename TImage::SizeType imageSize = img->GetLargestPossibleRegion().GetSize();

  std::cout << "Origin: ";
  for (typename TImage::SizeValueType j = 0; j < TImage::ImageDimension; j++)
    std::cout << milx::NumberToString(origin[j]) << "  ";
  std::cout << "\nSpacing: ";
  for (typename TImage::SizeValueType j = 0; j < TImage::ImageDimension; j++)
    std::cout << milx::NumberToString(spacing[j]) << "  ";
  std::cout << "\nSize: ";
  for (typename TImage::SizeValueType j = 0; j < TImage::ImageDimension; j++)
    std::cout << milx::NumberToString(imageSize[j]) << "  ";
  std::cout << "\nReal Size: ";
  for (typename TImage::SizeValueType j = 0; j < TImage::ImageDimension; j++)
    std::cout << milx::NumberToString(spacing[j] * imageSize[j]) << "  ";
  std::cout << "\nDirection/Orientation: " << std::endl;
  for (typename TImage::SizeValueType j = 0; j < TImage::ImageDimension; j++)
  {
    std::cout << "| ";
    for (typename TImage::SizeValueType k = 0; k < TImage::ImageDimension; k++)
    {
      std::cout << milx::NumberToString(direction(j, k));
      if (k < TImage::ImageDimension - 1)
        std::cout << ",\t";
    }
    std::cout << " |" << std::endl;
  }
}

template<class TImage>
double Image<TImage>::ImageMaximum(itk::SmartPointer<TImage> img)
{
  typedef itk::MinimumMaximumImageCalculator<TImage> ImageCalculatorFilterType;
  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(img);
  imageCalculatorFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  imageCalculatorFilter->Compute();

  return imageCalculatorFilter->GetMaximum();
}

template<class TImage>
double Image<TImage>::ImageMinimum(itk::SmartPointer<TImage> img)
{
  typedef itk::MinimumMaximumImageCalculator<TImage> ImageCalculatorFilterType;
  typename ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(img);
  imageCalculatorFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  imageCalculatorFilter->Compute();

  return imageCalculatorFilter->GetMinimum();
}

template<class TImage>
std::string Image<TImage>::ImageOrientation(itk::SmartPointer<TImage> img)
{
	std::map<itk::SpatialOrientation::ValidCoordinateOrientationFlags, std::string> codeToString;

	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL] = "AIL";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL] = "ASL";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI] = "RAI";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI] = "LAI";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS] = "RPS";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS] = "LPS";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP] = "RIP";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP] = "LIP";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP] = "RSP";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP] = "LSP";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA] = "RIA";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA] = "LIA";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA] = "RSA";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA] = "LSA";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP] = "IRP";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP] = "ILP";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP] = "SRP";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP] = "SLP";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA] = "IRA";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA] = "ILA";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA] = "SRA";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA] = "SLA";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI] = "RPI";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI] = "LPI";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS] = "RAS";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS] = "LAS";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI] = "PRI";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI] = "PLI";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI] = "ARI";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI] = "ALI";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS] = "PRS";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS] = "PLS";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS] = "ARS";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS] = "ALS";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR] = "IPR";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR] = "SPR";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR] = "IAR";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR] = "SAR";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL] = "IPL";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL] = "SPL";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL] = "IAL";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL] = "SAL";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR] = "PIR";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR] = "PSR";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR] = "AIR";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR] = "ASR";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL] = "PIL";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL] = "PSL";
	codeToString[itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_INVALID] = "Unknown";

	itk::SpatialOrientation::ValidCoordinateOrientationFlags orientFlag = itk::SpatialOrientationAdapter().FromDirectionCosines(img->GetDirection());
	std::string orientFlagStr = codeToString[orientFlag];

	return orientFlagStr;
}

//Filters
template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::CheckerBoard(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToCheckerBoard, const int numberOfSquares)
{
  typedef itk::CheckerBoardImageFilter<TImage> CheckerBoardType;
  typename CheckerBoardType::Pointer checkerBoardFilter = CheckerBoardType::New();
  checkerBoardFilter->SetInput1(img);
  checkerBoardFilter->SetInput2(imgToCheckerBoard);
  if(numberOfSquares != 0)
  {
    typename CheckerBoardType::PatternArrayType squares = checkerBoardFilter->GetCheckerPattern();
    squares.Fill(numberOfSquares);
    checkerBoardFilter->SetCheckerPattern(squares);
  }
  checkerBoardFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    checkerBoardFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Generating Checker Board");
    PrintError(ex.GetDescription());
  }

  return checkerBoardFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::RescaleIntensity(itk::SmartPointer<TImage> img, float minValue, float maxValue)
{
  typedef itk::RescaleIntensityImageFilter< TImage, TImage > RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(img);
  rescaleFilter->SetOutputMinimum(minValue);
  rescaleFilter->SetOutputMaximum(maxValue);
  rescaleFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    rescaleFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Rescaling Intensities");
    PrintError(ex.GetDescription());
  }

  return rescaleFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::InvertIntensity(itk::SmartPointer<TImage> img, float maxValue)
{
  typedef itk::InvertIntensityImageFilter<TImage, TImage> InvertIntensityImageFilterType;
  typename InvertIntensityImageFilterType::Pointer invert = InvertIntensityImageFilterType::New();
  invert->SetInput(img);
  invert->SetMaximum(maxValue);
  invert->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    invert->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Invert Intensity");
    PrintError(ex.GetDescription());
  }

  return invert->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::HistogramEqualisation(itk::SmartPointer<TImage> img, float alpha, float beta, float radius)
{
  typename TImage::SizeType radii;
  radii.Fill(radius);
  typedef itk::AdaptiveHistogramEqualizationImageFilter<TImage> EqualiseFilterType;
  typename EqualiseFilterType::Pointer equaliseFilter = EqualiseFilterType::New();
  equaliseFilter->SetInput(img);
  equaliseFilter->SetAlpha(alpha);
  equaliseFilter->SetBeta(beta);
  equaliseFilter->SetRadius(radii);
  equaliseFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    equaliseFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Histogram Equalisation");
    PrintError(ex.GetDescription());
  }

  return equaliseFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::GradientMagnitude(itk::SmartPointer<TImage> img)
{
  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<TImage, TImage> GradImageFilterType;
  typename GradImageFilterType::Pointer gradientMagnitudeFilter = GradImageFilterType::New();
  gradientMagnitudeFilter->SetInput(img);
  gradientMagnitudeFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    gradientMagnitudeFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Gradient Magnitude");
    PrintError(ex.GetDescription());
  }

  return gradientMagnitudeFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::SobelEdges(itk::SmartPointer<TImage> img)
{
  typedef itk::SobelEdgeDetectionImageFilter<TImage, TImage> SobelImageFilterType;
  typename SobelImageFilterType::Pointer sobelEdgeFilter = SobelImageFilterType::New();
  sobelEdgeFilter->SetInput(img);
  sobelEdgeFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    sobelEdgeFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Sobel Edges");
    PrintError(ex.GetDescription());
  }

  return sobelEdgeFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::Laplacian(itk::SmartPointer<TImage> img)
{
  typedef itk::LaplacianRecursiveGaussianImageFilter<TImage, TImage> LaplacianImageFilterType;
  typename LaplacianImageFilterType::Pointer laplacianEdgeFilter = LaplacianImageFilterType::New();
  laplacianEdgeFilter->SetInput(img);
  laplacianEdgeFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    laplacianEdgeFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Laplacian");
    PrintError(ex.GetDescription());
  }

  return laplacianEdgeFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::CannyEdges(itk::SmartPointer<TImage> img, float variance, float lowerThreshold, float upperThreshold)
{
  typedef itk::CannyEdgeDetectionImageFilter<TImage, TImage> CannyImageFilterType;
  typename CannyImageFilterType::Pointer cannyEdgeFilter = CannyImageFilterType::New();
  cannyEdgeFilter->SetInput(img);
  cannyEdgeFilter->SetVariance( variance );
  cannyEdgeFilter->SetUpperThreshold( upperThreshold );
  cannyEdgeFilter->SetLowerThreshold( lowerThreshold );
  cannyEdgeFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    cannyEdgeFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Canny Edges");
    PrintError(ex.GetDescription());
  }

  return cannyEdgeFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer< itk::Image<float, TImage::ImageDimension> > Image<TImage>::Normalization(itk::SmartPointer<TImage> img)
{
  typedef itk::NormalizeImageFilter< TImage, itk::Image<float, TImage::ImageDimension> > NormalizeImageFilterType;
  typename NormalizeImageFilterType::Pointer normalizeFilter = NormalizeImageFilterType::New();
  normalizeFilter->SetInput(img);
  normalizeFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    normalizeFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Normalization");
    PrintError(ex.GetDescription());
  }

  return normalizeFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::MatchInformation(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToMatch, bool originOnly)
{
  typedef itk::ChangeInformationImageFilter<TImage> ChangeInfoImageFilterType;
  typename ChangeInfoImageFilterType::Pointer change = ChangeInfoImageFilterType::New();
  change->SetInput(img);
  change->SetReferenceImage(imgToMatch);
  change->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  change->UseReferenceImageOn();
  if(!originOnly)
  {
    change->ChangeDirectionOn();
    change->ChangeSpacingOn();
    change->ChangeRegionOn();
  }
  else
  {
    change->ChangeDirectionOff();
    change->ChangeSpacingOff();
    change->ChangeRegionOff();
  }
  change->ChangeOriginOn();
  try
  {
    change->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Matching Info");
    PrintError(ex.GetDescription());
  }

  return change->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::MatchHistogram(itk::SmartPointer<TImage> img, itk::SmartPointer<TImage> imgToMatch, const int bins)
{
  typedef itk::HistogramMatchingImageFilter <TImage, TImage> MatchingFilterType;
  typename MatchingFilterType::Pointer matcher = MatchingFilterType::New();
  matcher->SetInput(img);
  matcher->SetReferenceImage(imgToMatch);
  matcher->SetNumberOfHistogramLevels(bins);
  matcher->SetNumberOfMatchPoints(7);
  matcher->ThresholdAtMeanIntensityOn();
  matcher->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    matcher->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Matching Histogram");
    PrintError(ex.GetDescription());
  }

  return matcher->GetOutput();
}

template<class TImage>
template<typename TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::DistanceMap(itk::SmartPointer<TImage> img, const bool binaryImage, const bool signedDistances, const bool computeInsideObject, const bool squaredDistance)
{
  typedef itk::ImageToImageFilter<TImage, TOutImage> ImageFilterType;

  typename ImageFilterType::Pointer filter;
  if(computeInsideObject)
  {
    typedef itk::ApproximateSignedDistanceMapImageFilter<TImage, TOutImage> DistanceMapImageFilterType;
    typename DistanceMapImageFilterType::Pointer distanceMapFilter = DistanceMapImageFilterType::New();
    filter = distanceMapFilter;
  }
  else if(signedDistances)
  {
    typedef itk::SignedMaurerDistanceMapImageFilter<TImage, TOutImage> DistanceMapImageFilterType;
    typename DistanceMapImageFilterType::Pointer distanceMapFilter = DistanceMapImageFilterType::New();
    distanceMapFilter->UseImageSpacingOn();
    if(!squaredDistance)
      distanceMapFilter->SquaredDistanceOff();
    else
      distanceMapFilter->SquaredDistanceOn();
    filter = distanceMapFilter;
  }
  else
  {
    typedef itk::DanielssonDistanceMapImageFilter<TImage, TOutImage> DistanceMapImageFilterType;
    typename DistanceMapImageFilterType::Pointer distanceMapFilter = DistanceMapImageFilterType::New();
    distanceMapFilter->UseImageSpacingOn();
    if(!binaryImage)
      distanceMapFilter->InputIsBinaryOff();
    else
      distanceMapFilter->InputIsBinaryOn();
    if(!squaredDistance)
      distanceMapFilter->SquaredDistanceOff();
    else
      distanceMapFilter->SquaredDistanceOn();
    filter = distanceMapFilter;
  }

//    std::cout << filter->GetNameOfClass() << std::endl;
  filter->SetInput(img);
  filter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    filter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Distance Map");
    PrintError(ex.GetDescription());
  }

  return filter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::FlipImage(itk::SmartPointer<TImage> img, bool xAxis, bool yAxis, bool zAxis, bool aboutOrigin)
{
  itk::FixedArray<bool, imgDimension> flipAxes;
  flipAxes[0] = xAxis;
  flipAxes[1] = yAxis;
  flipAxes[2] = zAxis;
//    flipAxes[3] = false;

  typedef itk::FlipImageFilter<TImage> FlipImageFilterType;
  typename FlipImageFilterType::Pointer flipFilter = FlipImageFilterType::New();
  flipFilter->SetInput(img);
  flipFilter->SetFlipAxes(flipAxes);
  if(aboutOrigin)
    flipFilter->FlipAboutOriginOn();
  else
    flipFilter->FlipAboutOriginOff();
  flipFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    flipFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Flip Image");
    PrintError(ex.GetDescription());
  }

  return flipFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::PadImageByConstant(itk::SmartPointer<TImage> img, size_t xAxis, size_t yAxis, size_t zAxis, typename TImage::PixelType value)
{
  typename TImage::SizeType lowerExtendRegion;
  lowerExtendRegion[0] = xAxis;
  lowerExtendRegion[1] = yAxis;
  lowerExtendRegion[2] = zAxis;

  typename TImage::SizeType upperExtendRegion;
  upperExtendRegion[0] = xAxis;
  upperExtendRegion[1] = yAxis;
  upperExtendRegion[2] = zAxis;

  typename TImage::PixelType constantPixel = value;

  typedef itk::ConstantPadImageFilter<TImage, TImage> ConstantPadImageFilterType;
  typename ConstantPadImageFilterType::Pointer padFilter = ConstantPadImageFilterType::New();
  padFilter->SetInput(img);
  padFilter->SetPadLowerBound(lowerExtendRegion);
  padFilter->SetPadUpperBound(upperExtendRegion);
  padFilter->SetConstant(constantPixel);
  padFilter->Update();
  padFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    padFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Padding Image");
    PrintError(ex.GetDescription());
  }

  return padFilter->GetOutput();
}

template<class TImage>
template<typename TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::AnisotropicDiffusion(itk::SmartPointer<TImage> img, const int iterations, const float timestep)
{
  typedef itk::GradientAnisotropicDiffusionImageFilter<TImage, TOutImage> GradientAnisotropicDiffusionImageFilterType;
  typename GradientAnisotropicDiffusionImageFilterType::Pointer gradAniDiff = GradientAnisotropicDiffusionImageFilterType::New();
  gradAniDiff->SetInput(img);
  gradAniDiff->SetUseImageSpacing(true);
  gradAniDiff->SetTimeStep(timestep);
  PrintInfo("Time step " + milx::NumberToString(gradAniDiff->GetTimeStep()));
  gradAniDiff->SetNumberOfIterations(iterations);
  gradAniDiff->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  PrintInfo("Conductance " + milx::NumberToString(gradAniDiff->GetConductanceParameter()));
  try
  {
    gradAniDiff->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Invert Intensity");
    PrintError(ex.GetDescription());
  }

  return gradAniDiff->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::Bilateral(itk::SmartPointer<TImage> img, const float sigmaRange, const float sigmaSpatial)
{
  typedef itk::BilateralImageFilter<TImage, TImage> BilateralImageFilterType;
  typename BilateralImageFilterType::Pointer bilateralFilter = BilateralImageFilterType::New();
  bilateralFilter->SetInput(img);
  bilateralFilter->SetRangeSigma(sigmaRange);
  bilateralFilter->SetDomainSigma(sigmaSpatial);
  bilateralFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    bilateralFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Bilateral Smoothing");
    PrintError(ex.GetDescription());
  }

  return bilateralFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::GaussianSmooth(itk::SmartPointer<TImage> img, const float variance)
{
  typedef itk::SmoothingRecursiveGaussianImageFilter<TImage, TImage> GuassianImageFilterType;
  typename GuassianImageFilterType::Pointer smoothGaussian = GuassianImageFilterType::New();
  smoothGaussian->SetInput(img);
  smoothGaussian->SetSigma(variance);
  smoothGaussian->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    smoothGaussian->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Gaussian Smoothing");
    PrintError(ex.GetDescription());
  }

  return smoothGaussian->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::Median(itk::SmartPointer<TImage> img, const int radius)
{
  typedef itk::MedianImageFilter<TImage, TImage> MedianImageFilterType;
  typename MedianImageFilterType::Pointer medianFilter = MedianImageFilterType::New();
  medianFilter->SetInput(img);
  typename MedianImageFilterType::InputSizeType radiusInput;
  radiusInput.Fill(radius);
  medianFilter->SetRadius(radiusInput);
  medianFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    medianFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Median Image");
    PrintError(ex.GetDescription());
  }

  return medianFilter->GetOutput();
}

//Labelling
template<class TImage>
template<typename TMaskImage>
itk::SmartPointer<TImage> Image<TImage>::MaskImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TMaskImage> maskImg)
{
  typedef itk::MaskImageFilter<TImage, TMaskImage> MaskFilterType;
  typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  maskFilter->SetInput(img);
  maskFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
#if ITK_MAJOR_VERSION > 3
  maskFilter->SetCoordinateTolerance(CoordinateTolerance);
  maskFilter->SetDirectionTolerance(DirectionTolerance);
  maskFilter->SetMaskImage(maskImg);
#else
  maskFilter->SetInput2(maskImg);
#endif

  try
  {
    maskFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Masking Image");
    PrintError(ex.GetDescription());
  }

  return maskFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::RelabelImage(itk::SmartPointer<TImage> labelledImg)
{
  typedef itk::ConnectedComponentImageFilter<TImage, TImage> ConnectedComponentImageFilterType;
  typedef itk::RelabelComponentImageFilter<TImage, TImage> RelabelImageType;

  typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
    connected->SetInput(labelledImg);

  typename RelabelImageType::Pointer relabelImageFilter = RelabelImageType::New();
    relabelImageFilter->SetInput(connected->GetOutput());

  try
  {
    relabelImageFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Relabelling Image");
    PrintError(ex.GetDescription());
  }

  return relabelImageFilter->GetOutput();
}

#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::BinaryContour(itk::SmartPointer<TImage> img, const float foregroundValue, const float backgroundValue)
{
  typedef itk::BinaryContourImageFilter<TImage, TImage> BinaryContourType;

  typename BinaryContourType::Pointer contourFilter = BinaryContourType::New();
  contourFilter->SetInput(img);
  contourFilter->SetForegroundValue(foregroundValue);
  contourFilter->SetBackgroundValue(backgroundValue);
  contourFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    contourFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Generating Contour");
    PrintError(ex.GetDescription());
  }

  return contourFilter->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::LabelContour(itk::SmartPointer<TImage> img, const bool fullyConnected, const float backgroundValue)
{
  typedef itk::LabelContourImageFilter<TImage, TImage> LabelContourType;

  typename LabelContourType::Pointer contourFilter = LabelContourType::New();
  contourFilter->SetInput(img);
  contourFilter->SetFullyConnected(fullyConnected);
  contourFilter->SetBackgroundValue(backgroundValue);
  contourFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    contourFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Generating Label Contours");
    PrintError(ex.GetDescription());
  }

  return contourFilter->GetOutput();
}

template<class TImage>
template<typename TMaskImage>
itk::SmartPointer<TImage> Image<TImage>::MaskAndCropImage(itk::SmartPointer<TImage> img, itk::SmartPointer<TMaskImage> maskImg, const size_t pixelPadding)
{
  typedef itk::LabelObject<unsigned char, TImage::ImageDimension> LabelObjectType;
  typedef itk::LabelMap<LabelObjectType> LabelMapType;

  typedef itk::LabelImageToLabelMapFilter<TMaskImage, LabelMapType> Image2LabelMapType;
  typename Image2LabelMapType::Pointer img2label = Image2LabelMapType::New();
    img2label->SetInput(maskImg);
    img2label->SetBackgroundValue(0);

  typedef itk::AutoCropLabelMapFilter<LabelMapType> AutoCropType;
  typename AutoCropType::Pointer autoCropper = AutoCropType::New();
    autoCropper->SetInput(img2label->GetOutput());
  typename AutoCropType::SizeType padSize = {{pixelPadding, pixelPadding, pixelPadding}};
    autoCropper->SetCropBorder(padSize);
    autoCropper->Update();

  typename AutoCropType::SizeType actualPadSize = autoCropper->GetCropBorder();
  milx::PrintDebug("Auto cropping with " + milx::NumberToString(actualPadSize[0]) + " pixel padding.");

  return ExtractSubImage<TImage>(img, autoCropper->GetRegion());
}

template<class TImage>
std::vector<unsigned char> Image<TImage>::LabelValues(itk::SmartPointer<TImage> labelledImg)
{
  typedef itk::LabelObject<unsigned char, milx::imgDimension> LabelObjectType;
  typedef itk::LabelMap<LabelObjectType>   LabelMapType;
  typedef itk::LabelImageToLabelMapFilter<TImage, LabelMapType> LabelImageToLabelMapType;

  typename LabelImageToLabelMapType::Pointer labelImageToLabelMapFilter = LabelImageToLabelMapType::New();
    labelImageToLabelMapFilter->SetInput(labelledImg);
    labelImageToLabelMapFilter->Update();
  typename LabelMapType::Pointer labelMap = labelImageToLabelMapFilter->GetOutput();

  std::vector<unsigned char> labelValues = labelMap->GetLabels();
  std::cout << "Found " << labelMap->GetNumberOfLabelObjects() << std::endl;

  return labelValues;
}

template<class TImage>
template<typename TLabelImage>
itk::SmartPointer< itk::Image<itk::RGBPixel<unsigned char>, TImage::ImageDimension> > Image<TImage>::Overlay(itk::SmartPointer<TImage> img, itk::SmartPointer<TLabelImage> overlayImg)
{
  itk::SmartPointer<TImage> rescaledImg = RescaleIntensity(img, 0.0, 255.0); ///Rescale intensities since images cast to 8-bit image

  typedef unsigned char charPixelType;
  typedef itk::Image<charPixelType, TImage::ImageDimension> charImageType;

  typedef itk::CastImageFilter<TLabelImage, charImageType> CastImageType;
  typename CastImageType::Pointer castFilter = CastImageType::New();
  castFilter->SetInput(overlayImg);
  castFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    castFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Casting in Overlay");
    PrintError(ex.GetDescription());
  }

  ///Cast Labels to unsigned char, since LabelOverlayImageFilter needs labels as integers
  typedef itk::LabelOverlayImageFilter< TImage, charImageType, itk::Image<itk::RGBPixel<unsigned char>, TImage::ImageDimension> > OverlayType;
  typename OverlayType::Pointer overlayFilter = OverlayType::New();
  overlayFilter->SetInput( rescaledImg );
  overlayFilter->SetLabelImage( castFilter->GetOutput() );
  overlayFilter->SetOpacity( 0.6 );
  overlayFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    overlayFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Generating Overlay");
    PrintError(ex.GetDescription());
  }

  return overlayFilter->GetOutput();
}

template<class TImage>
template<typename TLabelImage>
itk::SmartPointer< itk::Image<itk::RGBPixel<unsigned char>, TImage::ImageDimension> > Image<TImage>::OverlayContour(itk::SmartPointer<TImage> img, itk::SmartPointer<TLabelImage> overlayImg)
{
  itk::SmartPointer<TImage> rescaledImg = RescaleIntensity(img, 0.0, 255.0); ///Rescale intensities since images cast to 8-bit image

  typedef unsigned char charPixelType;
  typedef itk::Image<charPixelType, TImage::ImageDimension> charImageType;

  typedef itk::LabelContourImageFilter<TLabelImage, charImageType> LabelContourType;
  typename LabelContourType::Pointer contourFilter = LabelContourType::New();
  contourFilter->SetInput(overlayImg);
  contourFilter->SetFullyConnected(true);
  contourFilter->SetBackgroundValue(0.0);
  contourFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    contourFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Generating Label Contours in Overlay");
    PrintError(ex.GetDescription());
  }

  ///Cast Labels to unsigned char, since LabelOverlayImageFilter needs labels as integers
  typedef itk::LabelOverlayImageFilter< TImage, charImageType, itk::Image<itk::RGBPixel<unsigned char>, TImage::ImageDimension> > OverlayType;
  typename OverlayType::Pointer overlayFilter = OverlayType::New();
  overlayFilter->SetInput( rescaledImg );
  overlayFilter->SetLabelImage( contourFilter->GetOutput() );
  overlayFilter->SetOpacity( 0.6 );
  overlayFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    overlayFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Generating Overlay Contour");
    PrintError(ex.GetDescription());
  }

  return overlayFilter->GetOutput();
}
#endif // (ITK_REVIEW || ITK_VERSION_MAJOR > 3)

//Thresholding
template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::ThresholdAboveImage(itk::SmartPointer<TImage> img, float outsideValue, float aboveValue)
{
  typedef itk::ThresholdImageFilter<TImage> ThresholdImageFilterType;
  typename ThresholdImageFilterType::Pointer thresholding = ThresholdImageFilterType::New();
  thresholding->SetInput(img);
  thresholding->SetOutsideValue(outsideValue);
  thresholding->ThresholdAbove(aboveValue);
  thresholding->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    thresholding->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Threshold");
    PrintError(ex.GetDescription());
  }

  return thresholding->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::ThresholdBelowImage(itk::SmartPointer<TImage> img, float outsideValue, float belowValue)
{
  typedef itk::ThresholdImageFilter<TImage> ThresholdImageFilterType;
  typename ThresholdImageFilterType::Pointer thresholding = ThresholdImageFilterType::New();
  thresholding->SetInput(img);
  thresholding->SetOutsideValue(outsideValue);
  thresholding->ThresholdBelow(belowValue);
  thresholding->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    thresholding->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Threshold");
    PrintError(ex.GetDescription());
  }

  return thresholding->GetOutput();
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::ThresholdImage(itk::SmartPointer<TImage> img, float outsideValue, float belowValue, float aboveValue)
{
  typedef itk::ThresholdImageFilter<TImage> ThresholdImageFilterType;
  typename ThresholdImageFilterType::Pointer thresholding = ThresholdImageFilterType::New();
  thresholding->SetInput(img);
  thresholding->SetOutsideValue(outsideValue);
  thresholding->ThresholdOutside(belowValue, aboveValue);
  thresholding->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    thresholding->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Threshold");
    PrintError(ex.GetDescription());
  }

  return thresholding->GetOutput();
}

template<class TImage>
template<typename TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::BinaryThresholdImage(itk::SmartPointer<TImage> img, float outsideValue, float insideValue, float belowValue, float aboveValue)
{
  typedef itk::BinaryThresholdImageFilter<TImage, TOutImage> ThresholdImageFilterType;
  typename ThresholdImageFilterType::Pointer thresholding = ThresholdImageFilterType::New();
  thresholding->SetInput(img);
  thresholding->SetOutsideValue(outsideValue);
  thresholding->SetInsideValue(insideValue);
  thresholding->SetLowerThreshold(belowValue);
  thresholding->SetUpperThreshold(aboveValue);
  thresholding->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    thresholding->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Threshold");
    PrintError(ex.GetDescription());
  }

  return thresholding->GetOutput();
}

template<class TImage>
double Image<TImage>::OtsuThreshold(itk::SmartPointer<TImage> img, const int bins)
{
  typedef itk::OtsuThresholdImageFilter<TImage, TImage> OtsuThresholdImageCalculatorType;
  typename OtsuThresholdImageCalculatorType::Pointer otsuThresholdImageCalculator = OtsuThresholdImageCalculatorType::New();
  otsuThresholdImageCalculator->SetInput(img);
  otsuThresholdImageCalculator->SetNumberOfHistogramBins(bins);
  otsuThresholdImageCalculator->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    otsuThresholdImageCalculator->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Otsu Threshold");
    PrintError(ex.GetDescription());
  }

  //create binary image
  return otsuThresholdImageCalculator->GetThreshold();
}

template<class TImage>
template<typename TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::OtsuThresholdImage(itk::SmartPointer<TImage> img, const int bins)
{
  typedef itk::OtsuThresholdImageFilter<TImage, TOutImage> OtsuThresholdImageCalculatorType;
  typename OtsuThresholdImageCalculatorType::Pointer otsuThresholdImageCalculator = OtsuThresholdImageCalculatorType::New();
  otsuThresholdImageCalculator->SetInput(img);
  otsuThresholdImageCalculator->SetNumberOfHistogramBins(bins);
  otsuThresholdImageCalculator->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    otsuThresholdImageCalculator->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Otsu Threshold of Image");
    PrintError(ex.GetDescription());
  }
  PrintDebug("Otsu Threshold - " + milx::NumberToString(otsuThresholdImageCalculator->GetThreshold()));

  //create binary image
  return otsuThresholdImageCalculator->GetOutput();
}

template<class TImage>
template<typename TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::OtsuMultipleThresholdImage(itk::SmartPointer<TImage> img, const int bins, const int noOfLabels)
{
  typedef itk::OtsuMultipleThresholdsImageFilter<TImage, TOutImage> OtsuThresholdImageType;
  typename OtsuThresholdImageType::Pointer otsuThresholdImageCalculator = OtsuThresholdImageType::New();
  otsuThresholdImageCalculator->SetInput(img);
  otsuThresholdImageCalculator->SetNumberOfHistogramBins(bins);
  otsuThresholdImageCalculator->SetNumberOfThresholds(noOfLabels);
  otsuThresholdImageCalculator->SetLabelOffset(0);
  otsuThresholdImageCalculator->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    otsuThresholdImageCalculator->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Otsu Multiple Threshold");
    PrintError(ex.GetDescription());
  }
  //print thresholds
  typename OtsuThresholdImageType::ThresholdVectorType thresholds = otsuThresholdImageCalculator->GetThresholds();
  std::string thresholdStrs;
  for(size_t j = 0; j < thresholds.size(); j ++)
    thresholdStrs += milx::NumberToString(thresholds[j]) + ", ";
  PrintDebug("Otsu Thresholds - " + thresholdStrs);

  return otsuThresholdImageCalculator->GetOutput();
}

template<class TImage>
template<typename TScalarImage>
itk::SmartPointer<TScalarImage> Image<TImage>::VectorMagnitude(itk::SmartPointer<TImage> img)
{
#if (ITK_VERSION_MAJOR > 3)
  typedef itk::VectorMagnitudeImageFilter<TImage, TScalarImage> MagnitudeImageFilterType;
#else
  typedef itk::GradientToMagnitudeImageFilter<TImage, TScalarImage> MagnitudeImageFilterType;
#endif
  typename MagnitudeImageFilterType::Pointer magnitudeFilter = MagnitudeImageFilterType::New();
  magnitudeFilter->SetInput(img);
  magnitudeFilter->AddObserver(itk::ProgressEvent(), ProgressUpdates);
  try
  {
    magnitudeFilter->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Computing Magnitude of Image");
    PrintError(ex.GetDescription());
  }

  return magnitudeFilter->GetOutput();
}

//Collection members
template<class TImage>
void Image<TImage>::InformationCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
  {
    Information(images[j]);
    std::cout << std::endl;
  }
}

template<class TImage>
void Image<TImage>::MatchHistogramCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TImage> refImage, const int bins)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = MatchHistogram(images[j], refImage, bins);
}

template<class TImage>
void Image<TImage>::CheckerboardCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TImage> refImage, const int squares)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = CheckerBoard(images[j], refImage, squares);
}

template<class TImage>
void Image<TImage>::RescaleIntensityCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float belowValue, float aboveValue)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = RescaleIntensity(images[j], belowValue, aboveValue);
}

template<class TImage>
void Image<TImage>::InvertIntensityCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
  {
    float maxValue = ImageMaximum(images[j]);
    images[j] = InvertIntensity(images[j], maxValue);
  }
}

template<class TImage>
void Image<TImage>::FlipCollection(std::vector< typename itk::SmartPointer<TImage> > &images, bool xAxis, bool yAxis, bool zAxis, bool aboutOrigin)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
  {
    images[j] = FlipImage(images[j], xAxis, yAxis, zAxis, aboutOrigin);
  }
}

template<class TImage>
template<typename TOutImage>
std::vector< typename itk::SmartPointer<TOutImage> > Image<TImage>::AnisotropicDiffusionCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, const int iterations, float timestep)
{
  const size_t n = images.size();

  std::vector< itk::SmartPointer<TOutImage> > collection;
  for(size_t j = 0; j < n; j ++)
  {
    typename TImage::SizeType imgSize = images[j]->GetLargestPossibleRegion().GetSize();

    if(timestep < 0.0) //auto set timestep
    {
      size_t maxDimension = 0;
      for(size_t k = 0; k < imgSize.GetSizeDimension(); k ++)
      {
        if(imgSize[k] > maxDimension)
          maxDimension = imgSize[k];
      }
      milx::PrintDebug("Timestep will be reciprocal of " + milx::NumberToString(maxDimension));
      timestep = 2.0/maxDimension;
    }

    itk::SmartPointer<TOutImage> newImage = AnisotropicDiffusion<TOutImage>(images[j], iterations, timestep);
    collection.push_back(newImage);
  }

  return collection;
}

template<class TImage>
void Image<TImage>::BilateralCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float sigmaRange, float sigmaSpatial)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = Bilateral(images[j], sigmaRange, sigmaSpatial);
}

template<class TImage>
template<typename TOutImage>
std::vector< typename itk::SmartPointer<TOutImage> > Image<TImage>::CastCollection(const std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = images.size();

  std::vector< itk::SmartPointer<TOutImage> > collection;
  for(size_t j = 0; j < n; j ++)
  {
    itk::SmartPointer<TOutImage> newImage = CastImage<TOutImage>(images[j]);
    collection.push_back(newImage);
  }

  return collection;
}

template<class TImage>
void Image<TImage>::MedianCollection(std::vector< typename itk::SmartPointer<TImage> > &images, const int radius)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = Median(images[j], radius);
}

template<class TImage>
void Image<TImage>::GradientMagnitudeCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = GradientMagnitude(images[j]);
}

template<class TImage>
void Image<TImage>::LaplacianCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = Laplacian(images[j]);
}

template<class TImage>
template<typename TOutImage>
std::vector< typename itk::SmartPointer<TOutImage> > Image<TImage>::DistanceMapCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, const bool binaryImage, const bool signedDistances, const bool computeInsideObject, const bool squaredDistance)
{
  const size_t n = images.size();

  std::vector< itk::SmartPointer<TOutImage> > collection;
  for(size_t j = 0; j < n; j ++)
  {
    itk::SmartPointer<TOutImage> newImage = DistanceMap<TOutImage>(images[j], binaryImage, signedDistances, computeInsideObject, squaredDistance);
    collection.push_back(newImage);
  }

  return collection;
}

//Thresholding
template<class TImage>
void Image<TImage>::ThresholdAboveCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float aboveValue)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = ThresholdAboveImage(images[j], outsideValue, aboveValue);
}

template<class TImage>
void Image<TImage>::ThresholdBelowCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float belowValue)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = ThresholdBelowImage(images[j], outsideValue, belowValue);
}

template<class TImage>
void Image<TImage>::ThresholdCollection(std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float belowValue, float aboveValue)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = ThresholdImage(images[j], outsideValue, belowValue, aboveValue);
}

template<class TImage>
template<typename TOutImage>
std::vector< typename itk::SmartPointer<TOutImage> > Image<TImage>::BinaryThresholdCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, float outsideValue, float insideValue, float belowValue, float aboveValue)
{
  const size_t n = images.size();

  std::vector< itk::SmartPointer<TOutImage> > collection;
  for(size_t j = 0; j < n; j ++)
  {
    itk::SmartPointer<TOutImage> newImage = BinaryThresholdImage<TOutImage>(images[j], outsideValue, insideValue, belowValue, aboveValue);
    collection.push_back(newImage);
  }

  return collection;
}

template<class TImage>
template<typename TOutImage>
std::vector< typename itk::SmartPointer<TOutImage> > Image<TImage>::OtsuMultipleThresholdCollection(const std::vector< typename itk::SmartPointer<TImage> > &images, const int bins, const int noOfLabels)
{
  const size_t n = images.size();

  std::vector< itk::SmartPointer<TOutImage> > collection;
  for(size_t j = 0; j < n; j ++)
  {
    itk::SmartPointer<TOutImage> newImage = OtsuMultipleThresholdImage<TOutImage>(images[j], bins, noOfLabels);
    collection.push_back(newImage);
  }

  return collection;
}

template<class TImage>
void Image<TImage>::SubsampleCollection(std::vector< typename itk::SmartPointer<TImage> > &images, unsigned factor)
{
  const size_t n = images.size();

  typename TImage::SizeType factors;
  factors[0] = factor;
  factors[1] = factor;
  factors[2] = factor;

  for (size_t j = 0; j < n; j++)
  {
    images[j] = SubsampleImage(images[j], factors);
  }
}

#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
template<class TImage>
template<class TMaskImage>
void Image<TImage>::MaskAndCropCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TMaskImage> maskImage, const size_t pixelPadding)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = MaskAndCropImage<TMaskImage>(images[j], maskImage, pixelPadding);
}
#endif // (ITK_REVIEW || ITK_VERSION_MAJOR > 3)

template<class TImage>
template<class TMaskImage>
void Image<TImage>::MaskCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TMaskImage> maskImage)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = MaskImage<TMaskImage>(images[j], maskImage);
}

template<class TImage>
void Image<TImage>::RelabelCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = RelabelImage(images[j]);
}

template<class TImage>
template<class TRefImage>
void Image<TImage>::ResampleCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TRefImage> refImage)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = ResampleImage<TRefImage>(images[j], refImage);
}

template<class TImage>
template<class TRefImage>
void Image<TImage>::ResampleLabelCollection(std::vector< typename itk::SmartPointer<TImage> > &images, itk::SmartPointer<TRefImage> refImage)
{
  const size_t n = images.size();

  for(size_t j = 0; j < n; j ++)
    images[j] = ResampleLabel<TRefImage>(images[j], refImage);
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::AddCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = images.size();

  itk::SmartPointer<TImage> result = DuplicateImage(images[0]);
  for(size_t j = 1; j < n; j ++)
    result = AddImages(result, images[j]);

  return result;
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::SubtractCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = images.size();

  itk::SmartPointer<TImage> result = DuplicateImage(images[0]);
  for(size_t j = 1; j < n; j ++)
    result = SubtractImages(result, images[j]);

  return result;
}

template<class TImage>
template<class TOutImage>
itk::SmartPointer<TOutImage> Image<TImage>::AverageCollection(std::vector< typename itk::SmartPointer<TImage> > &images)
{
  const size_t n = images.size();

  itk::SmartPointer<TImage> sumImg = AddCollection(images);
  itk::SmartPointer<TOutImage> averageImg = ScaleImage<TOutImage>(sumImg, 1.0/n);

  return averageImg;
}

template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::AverageVectorCollection(std::vector< typename itk::SmartPointer<TImage> > &images, int numberOfComponents)
{
  const size_t n = images.size();

  itk::SmartPointer<TImage> sumImg = AddCollection(images);
  itk::SmartPointer<TImage> averageImg = ScaleVectorImage(sumImg, 1.0/n, numberOfComponents);

  return averageImg;
}

template<class TImage>
template<class TJoinedImage>
itk::SmartPointer<TJoinedImage> Image<TImage>::JoinCollection(std::vector< typename itk::SmartPointer<TImage> > &images, const double origin, const double spacing)
{
  const size_t n = images.size();

  typedef itk::JoinSeriesImageFilter<TImage, TJoinedImage> JoinSeriesFilterType;
  typename JoinSeriesFilterType::Pointer joinSeries = JoinSeriesFilterType::New();
  joinSeries->SetOrigin(origin);
  joinSeries->SetSpacing(spacing);
  for(size_t j = 0; j < n; j ++)
    joinSeries->PushBackInput(images[j]);
  try
  {
    joinSeries->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Joining Images");
    PrintError(ex.GetDescription());
  }

  return joinSeries->GetOutput();
}

#if (ITK_REVIEW || ITK_VERSION_MAJOR > 3)
template<class TImage>
itk::SmartPointer<TImage> Image<TImage>::MergeLabelledImages(std::vector< typename itk::SmartPointer<TImage> > &images, const size_t mergeType)
{
  const size_t n = images.size();

  typedef itk::LabelObject<unsigned char, milx::imgDimension> LabelObjectType;
  typedef itk::LabelMap<LabelObjectType>   LabelMapType;
  typedef itk::LabelImageToLabelMapFilter<TImage, LabelMapType> LabelImageToLabelMapType;
  typedef itk::MergeLabelMapFilter<LabelMapType> MergerType;
  typename MergerType::Pointer merger = MergerType::New();
  if(mergeType == 0)
    merger->SetMethod(MergerType::KEEP); ///do its best to keep the label unchanged, but if a label is already used in a previous label map
  else if(mergeType == 1)
    merger->SetMethod(MergerType::AGGREGATE); ///If the same label is found several times in the label maps, the label objects with the same label are merged
  else if(mergeType == 2)
    merger->SetMethod(MergerType::PACK); ///relabel all the label objects by order of processing
  else if(mergeType == 3)
    merger->SetMethod(MergerType::STRICT); ///keeps the labels unchanged and raises an exception if the same label is found in several images

  for (size_t i = 0; i < n; i ++)
  {
    typename LabelImageToLabelMapType::Pointer labelImageToLabelMapFilter = LabelImageToLabelMapType::New();
      labelImageToLabelMapFilter->SetInput(images[i]);
      labelImageToLabelMapFilter->Update();

    std::cout << "Found " << labelImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " in the labelled image " << i << std::endl;

    merger->SetInput(i, labelImageToLabelMapFilter->GetOutput());
  }

  try
  {
    merger->Update();
  }
  catch (itk::ExceptionObject & ex )
  {
    PrintError("Failed Merging Labelled Images");
    PrintError(ex.GetDescription());
  }

  //Convert to image again
  typedef itk::LabelMapToLabelImageFilter<LabelMapType, TImage> LabelMapToLabelImageType;
  typename LabelMapToLabelImageType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageType::New();
    labelMapToLabelImageFilter->SetInput(merger->GetOutput());
    labelMapToLabelImageFilter->Update();

  return labelMapToLabelImageFilter->GetOutput();
}
#endif

} //end namespace milx

#endif //__MILXIMAGE_H
