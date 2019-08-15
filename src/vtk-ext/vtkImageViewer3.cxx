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
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkImageActor.h>
#include <vtkImageMapToWindowLevelColors.h>
#if(VTK_MAJOR_VERSION > 5)
  #include <vtkStreamingDemandDrivenPipeline.h>
  #include <vtkImageMapper3D.h>
  #include <vtkInformation.h>
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
      #if(VTK_MAJOR_VERSION > 7)
        this->IV->GetInputAlgorithm()->UpdateInformation();
		vtkInformation *requests = this->IV->GetInputAlgorithm()->GetOutputInformation(0);
		int *updateExtent = vtkStreamingDemandDrivenPipeline::GetWholeExtent(this->IV->GetInputInformation());
		requests->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), updateExtent, 6);
		//vtkInformation* outInfo = this->IV->GetInputAlgorithm()->GetOutputInformation(0);
		//outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), updateExtent);
        this->IV->GetInputAlgorithm()->Update();
	  #elif(VTK_MAJOR_VERSION > 5)
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
      std::cout << "Window: " << newWindow << ", Level: " << newLevel << std::endl;
      this->IV->Render();
    }

  vtkImageViewer3 *IV;
  double InitialWindow;
  double InitialLevel;
};

class vtkImageViewer3CursorCallback : public vtkCommand
{
public:
  static vtkImageViewer3CursorCallback *New() { return new vtkImageViewer3CursorCallback; }
  vtkImageViewer3CursorCallback() : vtkCommand()
  {
    dataPicker = vtkPointPicker::New();
  }

  ~vtkImageViewer3CursorCallback()
  {
    if (dataPicker)
      dataPicker->Delete();
  }

  void Execute(vtkObject *caller,
    unsigned long event,
    void *vtkNotUsed(callData))
  {
    if (this->IV->GetInput() == NULL)
    {
      return;
    }

    //Get the interactor
    vtkInteractorStyleImage *isi =
      static_cast<vtkInteractorStyleImage *>(caller);

    int *size = this->IV->GetRenderWindow()->GetSize();

    if (this->IV->GetCursorEnabled())
    {
      if (dataPicker->Pick(isi->GetInteractor()->GetEventPosition()[0],
        isi->GetInteractor()->GetEventPosition()[1],
        0,  // always zero.
        isi->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer()))
      {
        double picked[3];
        dataPicker->GetMapperPosition(picked);
        //isi->GetInteractor()->GetPicker()->GetPickPosition(picked);
      #if(VTK_MAJOR_VERSION > 5)
        this->IV->GetCursor()->SetCenter(picked[0], picked[1], picked[2]);
      #else
        this->IV->GetCursor()->SetModelBounds(this->IV->GetInput()->GetBounds());
        this->IV->GetCursor()->SetFocalPoint(picked[0], picked[1], picked[2]);
      #endif
        this->IV->GetCursor()->Update();
        this->IV->GetRenderWindow()->Modified();
        this->IV->UpdateDisplayExtent();
      }
    }

    this->IV->Render();
  }

  vtkImageViewer3 *IV;
  vtkPointPicker *dataPicker;
};

//ImageViewer
vtkStandardNewMacro(vtkImageViewer3);

vtkImageViewer3::vtkImageViewer3()
{
  NeurologicalView = true;
  CursorEnabled = false;
  cursor = NULL;
#if(VTK_MAJOR_VERSION < 6)
  cursorMapper = NULL;
#endif
  cursorActor = NULL;
}

vtkImageViewer3::~vtkImageViewer3()
{
  if (cursor)
  {
    cursor->Delete();
    cursor = NULL;
  }
#if(VTK_MAJOR_VERSION < 6)
  if (cursorMapper)
  {
    cursorMapper->Delete();
    cursorMapper = NULL;
  }
#endif
  if (cursorActor)
  {
    cursorActor->Delete();
    cursorActor = NULL;
  }
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
    std::cout << "Using Neurological (head-first) Orientation" << std::endl;
    switch (this->SliceOrientation)
    {
      case vtkImageViewer3::SLICE_ORIENTATION_XY: //Axial
        cam->SetFocalPoint(0,0,0);
        cam->SetPosition(0,0,-1); // -1 if medical ?
        cam->SetViewUp(0,1,0);
        break;

      case vtkImageViewer3::SLICE_ORIENTATION_XZ: //Coronal
        cam->SetFocalPoint(0,0,0);
        cam->SetPosition(0,1,0); // 1 if medical ?
        cam->SetViewUp(0,0,1);
        break;

      case vtkImageViewer3::SLICE_ORIENTATION_YZ: //Saggital
        cam->SetFocalPoint(0,0,0);
        cam->SetPosition(-1,0,0); // -1 if medical ?
        cam->SetViewUp(0,0,1);
        break;
    }
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
        vtkCommand::WindowLevelEvent, cbk); //WindowLevel
      this->InteractorStyle->AddObserver(
        vtkCommand::StartWindowLevelEvent, cbk); //WindowLevel
      this->InteractorStyle->AddObserver(
        vtkCommand::ResetWindowLevelEvent, cbk); //WindowLevel
      cbk->Delete();
      vtkImageViewer3CursorCallback *cbk2 = vtkImageViewer3CursorCallback::New();
      cbk2->IV = this;
      this->InteractorStyle->AddObserver(
        vtkCommand::MiddleButtonPressEvent, cbk2); //Cursor
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

void vtkImageViewer3::SetSliceOrientation(int orientation)
{
  Superclass::SetSliceOrientation(orientation);
#if(VTK_MAJOR_VERSION > 5)
  if(CursorEnabled)
    cursorActor->GetCursorAlgorithm()->SetReslicePlaneNormal(this->SliceOrientation);
#endif
}

/*void vtkImageViewer2::SetSlice(int slice)
{
  Superclass::SetSlice(slice);

}*/

void vtkImageViewer3::EnableCursor()
{
#if(VTK_MAJOR_VERSION > 5)
  if (!cursor)
    cursor = vtkResliceCursor::New();
  cursor->SetImage(this->GetInput());
  cursor->SetThickMode(0);
  cursor->SetThickness(1, 1, 1);
  if(!CursorEnabled)
    {
      //Set as center of image to initialise
      /*double bounds[6];
      this->GetInput()->GetBounds(bounds);
      cursor->SetCenter( (bounds[1]-bounds[0])/2, (bounds[3]-bounds[2])/2, (bounds[5]-bounds[4])/2);*/
      cursor->SetCenter(this->GetInput()->GetOrigin());
    }

  if (!cursorActor)
    cursorActor = vtkResliceCursorActor::New();
  cursorActor->GetCursorAlgorithm()->SetResliceCursor(cursor);
  cursorActor->GetCursorAlgorithm()->SetReslicePlaneNormal(this->SliceOrientation);
#else
  if (!cursor)
    cursor = vtkCursor3D::New();
  cursor->SetModelBounds(this->GetInput()->GetBounds());
  cursor->AllOn();
  cursor->OutlineOff();

  if (!cursorMapper)
    cursorMapper = vtkPolyDataMapper::New();
  cursorMapper->SetInputConnection(cursor->GetOutputPort());

  if (!cursorActor)
    cursorActor = vtkActor::New();
  cursorActor->GetProperty()->SetColor(1, 0, 0);
  cursorActor->SetMapper(cursorMapper);
#endif
  cursor->Update();

  this->Renderer->AddActor(cursorActor);
  this->Modified();
  this->GetRenderWindow()->Modified();
  this->UpdateDisplayExtent();
  this->Render();
  CursorEnabled = true;
}

void vtkImageViewer3::UpdateCursor()
{
  if(CursorEnabled)
    {
      cursor->Update();
  #if(VTK_MAJOR_VERSION > 5)
      cursorActor->GetCursorAlgorithm()->SetReslicePlaneNormal(this->SliceOrientation);
  #endif
    }
}

void vtkImageViewer3::DisableCursor()
{
#if(VTK_MAJOR_VERSION > 5)
  this->Renderer->RemoveActor(cursorActor);
#else
  if (cursor)
    cursor->AllOff();
#endif
  CursorEnabled = false;
}

double* vtkImageViewer3::GetCursorFocalPoint()
{
#if(VTK_MAJOR_VERSION > 5)
  return cursor->GetCenter();
#else
  return cursor->GetFocalPoint();
#endif
}

void vtkImageViewer3::SetCursorFocalPoint(double *point)
{
#if(VTK_MAJOR_VERSION > 5)
  cursor->SetCenter(point);
#else
  cursor->SetFocalPoint(point);
#endif
}
