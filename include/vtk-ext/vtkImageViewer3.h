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
#ifndef __vtkImageViewer3_h
#define __vtkImageViewer3_h

#include <vtkObjectFactory.h> //For image display from VTK5+
#include <vtkImageViewer2.h> //For image display from VTK5+
#if(VTK_MAJOR_VERSION > 5)
  #include <vtkRenderingImageModule.h> // For export macro
  #define VTK_EXT_EXPORT VTKRENDERINGIMAGE_EXPORT
#else
  #define VTK_EXT_EXPORT VTK_RENDERING_EXPORT
#endif
#include <vtkInteractorStyleImage2.h> //For interactor like scanner

class VTK_EXT_EXPORT vtkImageViewer3 : public vtkImageViewer2
{
public:
  static vtkImageViewer3 *New();
  vtkTypeMacro(vtkImageViewer3, vtkImageViewer2);
  void PrintSelf(ostream& os, vtkIndent indent);

  //~ inline void SetDirection(floatImageType::DirectionType direction)
  //~ { Direction = direction;  }

  //Whether to use head first (neurological standard) or feet first orientation
  vtkSetMacro(NeurologicalView, bool);
  vtkGetMacro(NeurologicalView, bool);
  vtkBooleanMacro(NeurologicalView, bool);

protected:
  vtkImageViewer3();
  virtual ~vtkImageViewer3() {}

  virtual void UpdateOrientation();
  virtual void InstallPipeline();

  //~ floatImageType::DirectionType Direction;

  bool NeurologicalView;

#if(VTK_MAJOR_VERSION > 5)
  friend class vtkImageViewer3Callback;
#endif

private:
	vtkImageViewer3(const vtkImageViewer3&) {}  // Not implemented.
	void operator=(const vtkImageViewer3&) {}  // Not implemented.
};

#endif
