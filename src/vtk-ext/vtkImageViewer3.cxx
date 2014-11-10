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
#include "vtkImageViewer3.h"

#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkCommand.h>
#include <vtkImageData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageActor.h>
#include <vtkImageMapToWindowLevelColors.h>
#if(VTK_MAJOR_VERSION > 5)
  #include "vtkStreamingDemandDrivenPipeline.h"
  #include "vtkImageMapper3D.h"
#endif

//----------------------------------------------------------------------------
class vtkImageViewer3Callback : public vtkCommand
{
public:
  static vtkImageViewer3Callback *New() { return new vtkImageViewer3Callback; }

  void Execute(vtkObject *caller,
               unsigned long event,
               void *vtkNotUsed(callData))
    {
      if (this->IV->GetInput() == NULL)
        {
        return;
        }

      // Reset

      if (event == vtkCommand::ResetWindowLevelEvent)
        {
      #if(VTK_MAJOR_VERSION > 5)
        this->IV->GetInputAlgorithm()->UpdateInformation();
        vtkStreamingDemandDrivenPipeline::SetUpdateExtent(
          this->IV->GetInputInformation(),
          vtkStreamingDemandDrivenPipeline::GetWholeExtent(
            this->IV->GetInputInformation()));
        this->IV->GetInputAlgorithm()->Update();
      #else
        this->IV->GetInput()->UpdateInformation();
        this->IV->GetInput()->SetUpdateExtent
          (this->IV->GetInput()->GetWholeExtent());
        this->IV->GetInput()->Update();
      #endif
        double *range = this->IV->GetInput()->GetScalarRange();
        this->IV->SetColorWindow(range[1] - range[0]);
        this->IV->SetColorLevel(0.5 * (range[1] + range[0]));
        this->IV->Render();
        return;
        }

      // Start

      if (event == vtkCommand::StartWindowLevelEvent)
        {
        this->InitialWindow = this->IV->GetColorWindow();
        this->InitialLevel = this->IV->GetColorLevel();
        return;
        }

      // Adjust the window level here

      vtkInteractorStyleImage *isi =
        static_cast<vtkInteractorStyleImage *>(caller);

      int *size = this->IV->GetRenderWindow()->GetSize();
      double window = this->InitialWindow;
      double level = this->InitialLevel;

      // Compute normalized delta

      double dx = 4.0 *
        (isi->GetWindowLevelCurrentPosition()[0] -
         isi->GetWindowLevelStartPosition()[0]) / size[0];
      double dy = 4.0 *
        (isi->GetWindowLevelStartPosition()[1] -
         isi->GetWindowLevelCurrentPosition()[1]) / size[1];

      // Scale by current values

      if (fabs(window) > 0.01)
        {
        dx = dx * window;
        }
      else
        {
        dx = dx * (window < 0 ? -0.01 : 0.01);
        }
      if (fabs(level) > 0.01)
        {
        dy = dy * level;
        }
      else
        {
        dy = dy * (level < 0 ? -0.01 : 0.01);
        }

      // Abs so that direction does not flip

      if (window < 0.0)
        {
        dx = -1*dx;
        }
      if (level < 0.0)
        {
        dy = -1*dy;
        }

      // Compute new window level

      double newWindow = dx + window;
      double newLevel;
      newLevel = level - dy;

      // Stay away from zero and really

      if (fabs(newWindow) < 0.01)
        {
        newWindow = 0.01*(newWindow < 0 ? -1 : 1);
        }
      if (fabs(newLevel) < 0.01)
        {
        newLevel = 0.01*(newLevel < 0 ? -1 : 1);
        }

      this->IV->SetColorWindow(newWindow);
      this->IV->SetColorLevel(newLevel);
      this->IV->Render();
    }

  vtkImageViewer3 *IV;
  double InitialWindow;
  double InitialLevel;
};

//ImageViewer
vtkStandardNewMacro(vtkImageViewer3);

vtkImageViewer3::vtkImageViewer3()
{
  NeurologicalView = true;
}

void vtkImageViewer3::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkImageViewer3::UpdateOrientation()
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;

  if(!cam)
    return;

  if(!NeurologicalView)
  {
//    return this->Superclass::UpdateOrientation();
      switch (this->SliceOrientation)
      {
        case vtkImageViewer3::SLICE_ORIENTATION_XY: //Axial
          cam->SetFocalPoint(0,0,0);
          cam->SetPosition(0,0,1); // -1 if medical ?
          cam->SetViewUp(0,1,0);
          break;

        case vtkImageViewer3::SLICE_ORIENTATION_XZ: //Coronal
          cam->SetFocalPoint(0,0,0);
          cam->SetPosition(0,-1,0); // 1 if medical ?
          cam->SetViewUp(0,0,1);
          break;

        case vtkImageViewer3::SLICE_ORIENTATION_YZ: //Saggital
          cam->SetFocalPoint(0,0,0);
          cam->SetPosition(-1,0,0); // -1 if medical ?
          cam->SetViewUp(0,0,1);
          break;
      }
  }
  else
  {
    // Set the camera position
    /*
    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    matrix->Identity();

    //only want rotations in the plane
    for(int j = 0; j < 3; j ++)
    {
        matrix->SetElement(j, 0, Direction(j,0));
        matrix->SetElement(j, 1, Direction(j,1));
        matrix->SetElement(j, 2, Direction(j,2));
    }*/

    std::cout << "Using Neurological (head-first) Orientation" << std::endl;
    switch (this->SliceOrientation)
    {
      case vtkImageViewer3::SLICE_ORIENTATION_XY: //Axial
        cam->SetFocalPoint(0,0,0);
        cam->SetPosition(0,0,-1); // -1 if medical ?
        cam->SetViewUp(0,1,0);
        //~ matrix->SetElement(1, 1, -1);
        break;

      case vtkImageViewer3::SLICE_ORIENTATION_XZ: //Coronal
        cam->SetFocalPoint(0,0,0);
        cam->SetPosition(0,1,0); // 1 if medical ?
        cam->SetViewUp(0,0,1);
        //~ matrix->SetElement(1, 1, -1);
        break;

      case vtkImageViewer3::SLICE_ORIENTATION_YZ: //Saggital
        cam->SetFocalPoint(0,0,0);
        cam->SetPosition(-1,0,0); // -1 if medical ?
        cam->SetViewUp(0,0,1);
        //~ matrix->SetElement(1, 1, -1);
        break;
    }
    /*
    //Transform camera
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
      transform->Identity();
    transform->Concatenate(matrix);

    cam->SetUserTransform(transform);
    this->Renderer->ResetCamera();*/
  }
}

//----------------------------------------------------------------------------
void vtkImageViewer3::InstallPipeline()
{
  if (this->RenderWindow && this->Renderer)
    {
    this->RenderWindow->AddRenderer(this->Renderer);
    }

  if (this->Interactor)
    {
    if (!this->InteractorStyle)
      {
      vtkInteractorStyleImage2 *scannerInteractor = vtkInteractorStyleImage2::New();
        scannerInteractor->SetViewer(this);
      this->InteractorStyle = scannerInteractor;
      vtkImageViewer3Callback *cbk = vtkImageViewer3Callback::New();
      cbk->IV = this;
      this->InteractorStyle->AddObserver(
        vtkCommand::WindowLevelEvent, cbk);
      this->InteractorStyle->AddObserver(
        vtkCommand::StartWindowLevelEvent, cbk);
      this->InteractorStyle->AddObserver(
        vtkCommand::ResetWindowLevelEvent, cbk);
      cbk->Delete();
      }

    this->Interactor->SetInteractorStyle(this->InteractorStyle);
    this->Interactor->SetRenderWindow(this->RenderWindow);
    }

  if (this->Renderer && this->ImageActor)
    {
    this->Renderer->AddViewProp(this->ImageActor);
    }

  if (this->ImageActor && this->WindowLevel)
    {
  #if(VTK_MAJOR_VERSION > 5)
    this->ImageActor->GetMapper()->SetInputConnection(
      this->WindowLevel->GetOutputPort());
  #else
    this->ImageActor->SetInput(this->WindowLevel->GetOutput());
  #endif
    }
}
