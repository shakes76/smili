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
  #include <vtkResliceCursor.h> //For image display from VTK5+
  #include <vtkResliceCursorActor.h> //For image display from VTK5+
  #include <vtkResliceCursorPolyDataAlgorithm.h>
#else
  #include <vtkCursor3D.h> //For image display from VTK5+
  #include <vtkPolyDataMapper.h> //For image display from VTK5+
  #include <vtkActor.h> //For image display from VTK5+
#endif
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

  //Whether to use head first (neurological standard) or feet first orientation
  vtkSetMacro(NeurologicalView, bool);
  vtkGetMacro(NeurologicalView, bool);
  vtkBooleanMacro(NeurologicalView, bool);

  //Same as those in vtkImageViewer2, but handles cursor also
  virtual void SetSliceOrientation(int orientation);

  //Same as those in vtkImageViewer2, but handles cursor also
  //virtual void SetSlice(int s);

  //Enable the image cursor/crosshairs
#if(VTK_MAJOR_VERSION > 5)
  inline void SetCursor(vtkResliceCursor* newCursor)
  {
    if(cursor)
      cursor->Delete();
    cursor = newCursor;
    cursorActor->GetCursorAlgorithm()->SetResliceCursor(cursor);
  }
  inline vtkResliceCursor* GetCursor()
  {
    return cursor;
  }
  inline void SetCursorActor(vtkResliceCursorActor* newCursorActor)
  {
    if(cursorActor)
      cursorActor->Delete();
    cursorActor = newCursorActor;
  }
  inline vtkProp3D* GetCursorActor()
  {
    return cursorActor;
  }
#else
  inline void SetCursor(vtkCursor3D* newCursor)
  {
    if (cursor)
      cursor->Delete();
    cursor = newCursor;
  }
  inline vtkCursor3D* GetCursor()
  {
    return cursor; 
  }
  inline void SetCursorActor(vtkActor* newCursorActor)
  {
    if (cursorActor)
      cursorActor->Delete();
    cursorActor = newCursorActor;
  }
  inline vtkProp3D* GetCursorActor()
  {
    return cursorActor;
  }
#endif
  void EnableCursor();
  void DisableCursor();
  void UpdateCursor();
  vtkGetMacro(CursorEnabled, bool);

  double* GetCursorFocalPoint();
  void SetCursorFocalPoint(double *point);

protected:
  vtkImageViewer3();
  virtual ~vtkImageViewer3();

  virtual void UpdateOrientation();
  virtual void InstallPipeline();

  bool NeurologicalView;
  bool CursorEnabled;

#if(VTK_MAJOR_VERSION > 5)
  vtkResliceCursor *cursor;
  vtkResliceCursorActor *cursorActor;
#else
  vtkCursor3D *cursor;
  vtkPolyDataMapper *cursorMapper;
  vtkActor *cursorActor;
#endif

#if(VTK_MAJOR_VERSION > 5)
  friend class vtkImageViewer3Callback;
  friend class vtkImageViewer3CursorCallback;
#endif

private:
	vtkImageViewer3(const vtkImageViewer3&) {}  // Not implemented.
	void operator=(const vtkImageViewer3&) {}  // Not implemented.
};

#endif
