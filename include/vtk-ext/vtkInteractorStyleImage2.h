/*=========================================================================

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
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
#ifndef __vtkInteractorStyleImage2_h
#define __vtkInteractorStyleImage2_h

#include <vtkObjectFactory.h> //For image display from VTK5+
#if(VTK_MAJOR_VERSION > 5)
  #include <vtkInteractionStyleModule.h> // For export macro
  #define VTK_EXT_STYLE_EXPORT VTKINTERACTIONSTYLE_EXPORT
#else
  #define VTK_EXT_STYLE_EXPORT VTK_RENDERING_EXPORT
#endif
#include <vtkInteractorStyleImage.h> //For image display from VTK5+
#include <vtkImageViewer2.h> //For image display from VTK5+

class VTK_EXT_STYLE_EXPORT vtkInteractorStyleImage2 : public vtkInteractorStyleImage 
{
public:
  static vtkInteractorStyleImage2 *New();
  vtkTypeMacro(vtkInteractorStyleImage2, vtkInteractorStyleImage);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void OnMouseWheelForward();
  virtual void OnMouseWheelBackward();

  //Overloaded the Slice style of the original image interactor
  //~ void SetInteractionModeToImageSlicing();

  inline vtkImageViewer2* GetViewer() 
  { return viewer;  }
  inline void SetViewer(vtkImageViewer2 *view)
  { viewer = view;  }

protected:
  vtkInteractorStyleImage2() {}
  ~vtkInteractorStyleImage2() {}
    
  vtkImageViewer2 * viewer;

private:
  vtkInteractorStyleImage2(const vtkInteractorStyleImage2&);  // Not implemented.
  void operator=(const vtkInteractorStyleImage2&);  // Not implemented.
};

#endif
