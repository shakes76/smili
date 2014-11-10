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
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

#include "milxFile.h"
#include "milxImage.h"
#include "vtkImageViewer3.h"

#include <sstream>

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
  
  ///Open image
  VisualizingImageType::Pointer image = VisualizingImageType::New();
  if(!milx::File::OpenImage<VisualizingImageType>(fileName, image))
  {
    milx::PrintError("Could not open file.");
    return EXIT_FAILURE;
  }
    
  ///Reorient
  vtkSmartPointer<vtkImageData> imageVTK = milx::Image<VisualizingImageType>::ConvertITKImageToVTKImage(image);
  vtkSmartPointer<vtkImageData> imageVTKReoriented = milx::Image<VisualizingImageType>::ApplyOrientationToVTKImage(imageVTK, image);
  
  // Visualize
  vtkSmartPointer<vtkImageViewer3> imageViewer = vtkSmartPointer<vtkImageViewer3>::New();
  imageViewer->SetInput( imageVTKReoriented );
  imageViewer->GetRenderWindow()->SetSize( 500, 500 );
  //~ imageViewer->SetSlice(bounds[5]/2); //show middle of volume
  imageViewer->GetRenderer()->ResetCamera();
 
  // Set up an interactor that does not respond to mouse events
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  imageViewer->SetupInteractor( renderWindowInteractor );
  imageViewer->Render();
  
  //~ imageViewer->SetSliceOrientationToYZ(); //Coronal
  imageViewer->SetSliceOrientationToXZ(); //Saggital
 
  // Start the event loop
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
