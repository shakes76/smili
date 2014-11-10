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
#include "vtkInteractorStyleImage2.h"

#include <vtkRenderWindowInteractor.h>

//ImageViewer
vtkStandardNewMacro(vtkInteractorStyleImage2);

void vtkInteractorStyleImage2::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkInteractorStyleImage2::OnMouseWheelForward() 
{
  if(this->Interactor->GetShiftKey())
  {
    this->Superclass::OnMouseWheelForward();
  }
  else
  {
    if(!viewer)
      return;

    if(viewer->GetSlice()+1 <= viewer->GetSliceMax())
    {
      viewer->SetSlice(viewer->GetSlice()+1);
    }

    this->StartState(VTKIS_NONE);
  }
}

void vtkInteractorStyleImage2::OnMouseWheelBackward()
{
  if(this->Interactor->GetShiftKey())
  {
    this->Superclass::OnMouseWheelBackward();
  }
  else
  {
    if(!viewer)
      return;

    if (viewer->GetSlice()-1 >= viewer->GetSliceMin())
    {
      viewer->SetSlice(viewer->GetSlice()-1);
    }

    this->StartState(VTKIS_NONE);
  }
}

//~ void vtkInteractorStyleImage2::SetInteractionModeToImageSlicing()
//~ {
  
//~ }
