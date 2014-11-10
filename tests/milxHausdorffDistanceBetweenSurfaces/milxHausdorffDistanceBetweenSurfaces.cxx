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
///Other Headers
#include <string>
#include <cmath>
#include <ctime>
///SMILI Headers
#include <milxFile.h>
#include <milxImage.h>
#include <milxModel.h>
///ITK
#include <itkLinearInterpolateImageFunction.h>
///VTK
#include <vtkFloatArray.h>

/**
 * Compute Hausdorff distances using distance maps and the surfaces.
 */
typedef float InputPixelType;
typedef itk::Image<InputPixelType, milx::imgDimension> InputImageType;

vtkSmartPointer<vtkFloatArray> SurfaceScalarsFromImage(vtkSmartPointer<vtkPolyData> surface, itk::SmartPointer<InputImageType> img, const bool absoluteValues = false);

int main(int argc, char *argv[])
{
    srand( time(NULL) );

    if (argc < 6)
    {
        cerr << "milxHausdorffDistanceBetweenSurfaces:" << endl;
        cerr << "Input Surface 1 Filename " << endl;
        cerr << "Input (Voxelised) Surface 1 Binary Image Filename " << endl;
        cerr << "Input Surface 2 Filename " << endl;
        cerr << "Input (Voxelised) Surface 2 Binary Image Filename " << endl;
        cerr << "Output Surface Filename " << endl;
        return EXIT_FAILURE;
    }

    std::string inputSurfaceFilename1 = argv[1];
    std::string inputImageFilename1 = argv[2];
    std::string inputSurfaceFilename2 = argv[3];
    std::string inputImageFilename2 = argv[4];
    std::string outputImageFilename = argv[5];

    ///Read Surface files
    vtkSmartPointer<vtkPolyData> surface1;
    vtkSmartPointer<vtkPolyData> surface2;

    milx::PrintInfo("Reading Surfaces...");
    if(!milx::File::OpenModel(inputSurfaceFilename1, surface1) || !milx::File::OpenModel(inputSurfaceFilename2, surface2))
    {
        std::cerr << "Error Reeading Surfaces." << std::endl;
        exit(EXIT_FAILURE);
    }

    ///Read images
    InputImageType::Pointer image1;
    InputImageType::Pointer image2;

    milx::PrintInfo("Reading Images...");
    if(!milx::File::OpenImage<InputImageType>(inputImageFilename1, image1) || !milx::File::OpenImage<InputImageType>(inputImageFilename2, image2))
    {
        std::cerr << "Error Reeading Images." << std::endl;
        exit(EXIT_FAILURE);
    }

    ///Compute Hausdorff distance using distance maps
    ///Compute DTs
    InputImageType::Pointer distanceImage1;
    InputImageType::Pointer distanceImage2;
    const bool binary = false, signedDistance = true, insideDistance = false, squaredDistance = false;

    milx::PrintInfo("Computing Distance maps...");
    distanceImage1 = milx::Image<InputImageType>::DistanceMap<InputImageType>(image1, binary, signedDistance, insideDistance, squaredDistance);
    distanceImage2 = milx::Image<InputImageType>::DistanceMap<InputImageType>(image2, binary, signedDistance, insideDistance, squaredDistance);

    ///Save Fwd and Bwd distance on meshes
//    vtkSmartPointer<vtkPolyData> hausdorffSurface1 = vtkSmartPointer<vtkPolyData>::New();
//    hausdorffSurface1->DeepCopy(surface1);
//    vtkSmartPointer<vtkPolyData> hausdorffSurface2 = vtkSmartPointer<vtkPolyData>::New();
//    hausdorffSurface2->DeepCopy(surface2);

    milx::PrintInfo("Computing Scalars for Surfaces...");
    const bool absValues = true;
    vtkSmartPointer<vtkFloatArray> scalars1 = SurfaceScalarsFromImage(surface1, distanceImage2, absValues);
    vtkSmartPointer<vtkFloatArray> scalars2 = SurfaceScalarsFromImage(surface2, distanceImage1, absValues);

    double range[2];
    scalars1->GetRange(range);
    std::cout << "Hausdorff Distance 1: " << range[1] << std::endl;
    scalars2->GetRange(range);
    std::cout << "Hausdorff Distance 2: " << range[1] << std::endl;

    surface1->GetPointData()->SetScalars(scalars1);
    surface2->GetPointData()->SetScalars(scalars2);

    ///Hausdorff as mean of fwd and bwd distances
    vtkSmartPointer<vtkPolyData> hausdorffSurface = surface1;

    ///Write results
    milx::PrintInfo("Writing Surfaces...");
    milx::File::SaveModel(outputImageFilename, hausdorffSurface);

    return EXIT_SUCCESS;
}

vtkSmartPointer<vtkFloatArray> SurfaceScalarsFromImage(vtkSmartPointer<vtkPolyData> surface, itk::SmartPointer<InputImageType> img, const bool absoluteValues)
{
  std::cerr << "Marking Surface" << std::endl;
  const int numberOfPoints = surface->GetNumberOfPoints();

  typedef itk::Point<double, 3> InputImagePointType;
  typedef itk::ContinuousIndex<double, 3 > ContinuousIndexType;

  typedef itk::LinearInterpolateImageFunction<InputImageType, double> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage(img);

  vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
  scalars->SetName("Distance");
  scalars->SetNumberOfTuples(numberOfPoints);
  scalars->SetNumberOfComponents(1);

  InputImagePointType point;
  ContinuousIndexType index;
  for(int i = 0; i < numberOfPoints; i++)
  {
    double position[3];
    surface->GetPoint(i, position);

    point[0] = position[0];
    point[1] = position[1];
    point[2] = position[2];

    img->TransformPhysicalPointToContinuousIndex(point, index);

    // If inside, mark
    if(interpolator->IsInsideBuffer(index))
    {

      double valueFound = 0.0;
      if(absoluteValues)
        valueFound = fabs(interpolator->EvaluateAtContinuousIndex(index));
      else
        valueFound = interpolator->EvaluateAtContinuousIndex(index);
      scalars->SetTuple1(i, valueFound);
    }
  }

  return scalars;
}
