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

#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkImageMapper.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageProperty.h>
#include <vtkImageSlice.h>

#include "milxFile.h"
#include "milxImage.h"
 
int main(int argc, char *argv[])
{
  if (argc < 2)
  {
      cerr << "View an image with correct orientation." << endl;
      cerr << "Usage: " << argv[0] << " <image>" << endl;
      exit(EXIT_FAILURE);
  }
  std::string fileName = argv[1];
  
  typedef itk::Image<float, 3> VisualizingImageType;
  
  //Control
  /*
    Left Mouse button triggers window level events
    SHIFT Left Mouse pans the camera
    CTRL SHIFT Left Mouse dollys (a positional zoom) the camera
    Middle mouse button pans the camera
    Right mouse button dollys the camera.
    SHIFT Right Mouse triggers pick events
  
    CTRL Left Mouse slices through the image
    SHIFT Middle Mouse slices through the image
    CTRL Right Mouse spins the camera
  
    R Reset the Window/Level
    X Reset to a sagittal view
    Y Reset to a coronal view
    Z Reset to an axial view
  */
  
  ///Open image
  VisualizingImageType::Pointer image = VisualizingImageType::New();
  if(!milx::File::OpenImage<VisualizingImageType>(fileName, image))
  {
    milx::PrintError("Could not open file.");
    return EXIT_FAILURE;
  }
  
  ///Reorient
  vtkSmartPointer<vtkImageData> imageVTK = milx::Image<VisualizingImageType>::ConvertITKImageToVTKImage(image);
  vtkSmartPointer<vtkImageData> imageVTKReoriented = milx::Image<VisualizingImageType>::ApplyOrientationToVTKImage(imageVTK, image, true);
 
  vtkSmartPointer<vtkImageResliceMapper> imageResliceMapper = vtkSmartPointer<vtkImageResliceMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  imageResliceMapper->SetInputConnection(imageVTKReoriented->GetProducerPort());
#else
  imageResliceMapper->SetInputData(imageVTKReoriented);
#endif
  imageResliceMapper->SliceFacesCameraOn();
  imageResliceMapper->SliceAtFocalPointOn();
  
  vtkSmartPointer<vtkImageProperty> ip = vtkSmartPointer<vtkImageProperty>::New();
    ip->SetAmbient(0.0);
    ip->SetDiffuse(1.0);
    ip->SetOpacity(1.0);
    ip->SetInterpolationTypeToLinear();
 
  vtkSmartPointer<vtkImageSlice> imageSlice = vtkSmartPointer<vtkImageSlice>::New();
    imageSlice->SetMapper(imageResliceMapper);
    imageSlice->SetProperty(ip);
 
  // Setup renderers
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddViewProp(imageSlice);
  renderer->ResetCamera();
 
  // Setup render window
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(512, 512);
  renderWindow->AddRenderer(renderer);
 
  // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
 
  vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    //~ style->SetInteractionModeToImage3D();
    style->SetInteractionModeToImageSlicing();
 
  renderWindowInteractor->SetInteractorStyle(style);
 
  // Render and start interaction
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Initialize();
 
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
}
