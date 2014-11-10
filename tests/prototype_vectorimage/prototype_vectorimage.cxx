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
#include <itkImage.h>
#include <itkImageFileReader.h>
#include "itkImageFileWriter.h"
#include <itkExtractImageFilter.h>
#include "itkComposeImageFilter.h"
#include "itkVectorImage.h"
 
int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
    std::cerr << "Convert a 4D image to a vector image." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile outputImageFile" << std::endl;
    return EXIT_FAILURE;
    }
 
  typedef float PixelType;
  typedef itk::Image<PixelType, 4>         ImageType;
  typedef itk::Image<PixelType, 3>         Image3DType;
  typedef itk::VectorImage<PixelType, 3>         VectorImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
 
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  ImageType::Pointer image = reader->GetOutput();
    
  ImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
  std::cout << "Size: " << size[0] << "x" << size[1] << "x" << size[2] << "x" << size[3] << std::endl;
    
  typedef itk::ComposeImageFilter<Image3DType> ImageToVectorImageFilterType;
  ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
  for(size_t j = 0; j < size[3]; j ++)
  {
    ImageType::IndexType desiredStart;
    desiredStart.Fill(0);
    desiredStart[3] = j;
   
    ImageType::SizeType desiredSize;
    desiredSize[0] = size[0];
    desiredSize[1] = size[1];
    desiredSize[2] = size[2];
    desiredSize[3] = 0;
   
    ImageType::RegionType desiredRegion(desiredStart, desiredSize);
    
    //Extract 3D image 'slice' of 4D volume
    typedef itk::ExtractImageFilter<ImageType, Image3DType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetExtractionRegion(desiredRegion);
    filter->SetInput(image);
  #if ITK_VERSION_MAJOR >= 4
    filter->SetDirectionCollapseToIdentity(); // This is required.
  #endif
    filter->Update();
      
    //add to vector image
    imageToVectorImageFilter->SetInput(j, filter->GetOutput());
  }
  imageToVectorImageFilter->Update();
  
  typedef  itk::ImageFileWriter<VectorImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(imageToVectorImageFilter->GetOutput());
  writer->Update();
 
  return EXIT_SUCCESS;
}
