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
///VTK
#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkImageData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCellPicker.h>
#include <vtkOutlineFilter.h>
#include <vtkProperty.h>
#include <vtkImageActor.h>
#include <vtkImageMapToColors.h>
#include <vtkImageOrthoPlanes.h>
#include <vtkImagePlaneWidget.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
//milx
#include "milxFile.h"
#include "milxImage.h"

//----------------------------------------------------------------------------
class vtkOrthoPlanesCallback : public vtkCommand
{
public:
  static vtkOrthoPlanesCallback *New()
  { return new vtkOrthoPlanesCallback; }

  void Execute( vtkObject *caller, unsigned long vtkNotUsed( event ),
                void *callData )
  {
    vtkImagePlaneWidget* self =
      reinterpret_cast< vtkImagePlaneWidget* >( caller );
    if(!self) return;

    double* wl = static_cast<double*>( callData );

    if ( self == this->WidgetX )
      {
      this->WidgetY->SetWindowLevel(wl[0],wl[1],1);
      this->WidgetZ->SetWindowLevel(wl[0],wl[1],1);
      }
    else if( self == this->WidgetY )
      {
      this->WidgetX->SetWindowLevel(wl[0],wl[1],1);
      this->WidgetZ->SetWindowLevel(wl[0],wl[1],1);
      }
    else if (self == this->WidgetZ)
      {
      this->WidgetX->SetWindowLevel(wl[0],wl[1],1);
      this->WidgetY->SetWindowLevel(wl[0],wl[1],1);
      }
  }

  vtkOrthoPlanesCallback():WidgetX( 0 ), WidgetY( 0 ), WidgetZ ( 0 ) {}

  vtkImagePlaneWidget* WidgetX;
  vtkImagePlaneWidget* WidgetY;
  vtkImagePlaneWidget* WidgetZ;
};

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "prototype_image_planes test:" << endl;
        cerr << "Image Filename " << endl;
        return EXIT_FAILURE;
    }

    std::string inputSurfaceFilename = argv[1];

    //Load the vertex table from CSV file
    cout << "Reading file: " << inputSurfaceFilename << endl;
    typedef itk::Image< float, 3 >         ImageType;
    ImageType::Pointer image = ImageType::New();
    milx::File::OpenImage<ImageType>(inputSurfaceFilename, image);
    vtkSmartPointer<vtkImageData> imageData = milx::Image<ImageType>::ConvertITKImageToVTKImage(image);

    vtkSmartPointer<vtkOutlineFilter> outline = vtkSmartPointer<vtkOutlineFilter>::New();
      outline->SetInputConnection(imageData->GetProducerPort());

    vtkSmartPointer<vtkPolyDataMapper> outlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      outlineMapper->SetInputConnection(outline->GetOutputPort());

    vtkSmartPointer<vtkActor> outlineActor = vtkSmartPointer<vtkActor>::New();
      outlineActor->SetMapper( outlineMapper);

    vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> ren2 = vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
      renWin->SetMultiSamples(0);
      renWin->AddRenderer(ren2);
      renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);

    vtkSmartPointer<vtkCellPicker> picker =
      vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.005);

    vtkSmartPointer<vtkProperty> ipwProp =
      vtkSmartPointer<vtkProperty>::New();
    //assign default props to the ipw's texture plane actor

    vtkSmartPointer<vtkImagePlaneWidget> planeWidgetX = vtkSmartPointer<vtkImagePlaneWidget>::New();
      planeWidgetX->SetInteractor( iren);
      planeWidgetX->SetKeyPressActivationValue('x');
      planeWidgetX->SetPicker(picker);
      planeWidgetX->RestrictPlaneToVolumeOn();
      planeWidgetX->GetPlaneProperty()->SetColor(1,0,0);
      planeWidgetX->SetTexturePlaneProperty(ipwProp);
      planeWidgetX->TextureInterpolateOff();
      planeWidgetX->SetResliceInterpolateToNearestNeighbour();
      planeWidgetX->SetInput(imageData);
      planeWidgetX->SetPlaneOrientationToXAxes();
      planeWidgetX->SetSliceIndex(32);
      planeWidgetX->DisplayTextOn();
      planeWidgetX->On();
      planeWidgetX->InteractionOff();
      planeWidgetX->InteractionOn();

    vtkSmartPointer<vtkImagePlaneWidget> planeWidgetY = vtkSmartPointer<vtkImagePlaneWidget>::New();
      planeWidgetY->SetInteractor( iren);
      planeWidgetY->SetKeyPressActivationValue('y');
      planeWidgetY->SetPicker(picker);
      planeWidgetY->GetPlaneProperty()->SetColor(1,1,0);
      planeWidgetY->SetTexturePlaneProperty(ipwProp);
      planeWidgetY->TextureInterpolateOn();
      planeWidgetY->SetResliceInterpolateToLinear();
      planeWidgetY->SetInput(imageData);
      planeWidgetY->SetPlaneOrientationToYAxes();
      planeWidgetY->SetSlicePosition(102.4);
      planeWidgetY->SetLookupTable( planeWidgetX->GetLookupTable());
      planeWidgetY->DisplayTextOff();
      planeWidgetY->UpdatePlacement();
      planeWidgetY->On();

    vtkSmartPointer<vtkImagePlaneWidget> planeWidgetZ = vtkSmartPointer<vtkImagePlaneWidget>::New();
      planeWidgetZ->SetInteractor( iren);
      planeWidgetZ->SetKeyPressActivationValue('z');
      planeWidgetZ->SetPicker(picker);
      planeWidgetZ->GetPlaneProperty()->SetColor(0,0,1);
      planeWidgetZ->SetTexturePlaneProperty(ipwProp);
      planeWidgetZ->TextureInterpolateOn();
      planeWidgetZ->SetResliceInterpolateToCubic();
      planeWidgetZ->SetInput(imageData);
      planeWidgetZ->SetPlaneOrientationToZAxes();
      planeWidgetZ->SetSliceIndex(25);
      planeWidgetZ->SetLookupTable( planeWidgetX->GetLookupTable());
      planeWidgetZ->DisplayTextOn();
      planeWidgetZ->On();

    vtkSmartPointer<vtkImageOrthoPlanes> orthoPlanes = vtkSmartPointer<vtkImageOrthoPlanes>::New();
      orthoPlanes->SetPlane(0, planeWidgetX);
      orthoPlanes->SetPlane(1, planeWidgetY);
      orthoPlanes->SetPlane(2, planeWidgetZ);
      orthoPlanes->ResetPlanes();

    vtkSmartPointer<vtkOrthoPlanesCallback> cbk = vtkSmartPointer<vtkOrthoPlanesCallback>::New();
      cbk->WidgetX = planeWidgetX;
      cbk->WidgetY = planeWidgetY;
      cbk->WidgetZ = planeWidgetZ;
      planeWidgetX->AddObserver( vtkCommand::EndWindowLevelEvent, cbk );
      planeWidgetY->AddObserver( vtkCommand::EndWindowLevelEvent, cbk );
      planeWidgetZ->AddObserver( vtkCommand::EndWindowLevelEvent, cbk );
      
    double wl[2];
    planeWidgetZ->GetWindowLevel(wl);

    // Add a 2D image to test the GetReslice method
    vtkSmartPointer<vtkImageMapToColors> colorMap = vtkSmartPointer<vtkImageMapToColors>::New();
      colorMap->PassAlphaToOutputOff();
      colorMap->SetActiveComponent(0);
      colorMap->SetOutputFormatToLuminance();
      colorMap->SetInput(planeWidgetZ->GetResliceOutput());
      colorMap->SetLookupTable(planeWidgetX->GetLookupTable());

    vtkSmartPointer<vtkImageActor> imageActor = vtkSmartPointer<vtkImageActor>::New();
      imageActor->PickableOff();
      imageActor->SetInput(colorMap->GetOutput());

    // Visualize
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderer->AddActor(imageActor);
    renderer->SetBackground(1,1,1); // Background color white

    renderWindow->Render();
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}
