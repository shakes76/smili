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
// .NAME vtkDiffusionTensorGlyphFilter - control the generation and placement of glyphs at input points

#ifndef __vtkDiffusionTensorGlyphFilter_h
#define __vtkDiffusionTensorGlyphFilter_h

#include <vector>

#include "vtkSmartPointer.h"
#include "vtkProgrammableGlyphFilter.h"
#if(VTK_MAJOR_VERSION > 5)
  #include <vtkFiltersProgrammableModule.h> // For export macro
  #define VTK_GRAPHICS_EXPORT VTKFILTERSPROGRAMMABLE_EXPORT
#endif
#include "vtkSphereSource.h"
#include "vtkPolyData.h"

#include <itkVectorImage.h>

class VTK_GRAPHICS_EXPORT vtkDiffusionTensorGlyphFilter : public vtkProgrammableGlyphFilter
{
public:
  vtkTypeMacro(vtkDiffusionTensorGlyphFilter,vtkProgrammableGlyphFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  typedef itk::VectorImage<float, 3> VectorImageType;

  // Description
  static vtkDiffusionTensorGlyphFilter *New();

  void SetTensorImage(itk::SmartPointer<VectorImageType> img)
  { m_TensorImage = img;  }
  itk::SmartPointer<VectorImageType> GetTensorImage()
  { return m_TensorImage;  }
  void SetTensorImageIndices(std::vector<VectorImageType::IndexType> indices)
  { m_Indices = indices;  }
  std::vector<VectorImageType::IndexType> GetTensorImageIndices()
  { return m_Indices;  }
  void SetNormaliseVector(VectorImageType::PixelType vec)
  { m_NormVector = vec;  }
  VectorImageType::PixelType GetNormaliseVector()
  { return m_NormVector;  }
  void SetResolution(const unsigned res)
  { m_Resolution = res;  }
  unsigned GetResolution()
  { return m_Resolution;  }

  static double computeAmplitude(std::vector<double> SH, double x, double y, double z, int lmax);
  static int NforL (int lmax) { return ((lmax+1)*(lmax+2)/2); }
  static int LforN (int N) { return (2*(((int) (sqrt((float) (1+8*N)))-3)/4)); }
  static int index (int l, int m) { return (l*(l+1)/2 + m); }

protected:
  vtkDiffusionTensorGlyphFilter();
  ~vtkDiffusionTensorGlyphFilter();

  ///Overide the base class method to customise for the diffusion glyphs, colouring etc.
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//  virtual int FillInputPortInformation(int, vtkInformation *);

  VectorImageType::Pointer m_TensorImage;
  std::vector<VectorImageType::IndexType> m_Indices;
  VectorImageType::PixelType m_NormVector;
  unsigned m_Resolution;

private:
  vtkDiffusionTensorGlyphFilter(const vtkDiffusionTensorGlyphFilter&);  // Not implemented.
  void operator=(const vtkDiffusionTensorGlyphFilter&);  // Not implemented.
};

#endif
