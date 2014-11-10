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
#ifndef __vtkAreaSimplificationMetric_h
#define __vtkAreaSimplificationMetric_h

#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkReebGraphSimplificationMetric.h>

class VTK_RENDERING_EXPORT vtkAreaSimplificationMetric : public vtkReebGraphSimplificationMetric
{

public:
  vtkTypeMacro(vtkAreaSimplificationMetric, vtkReebGraphSimplificationMetric);

  static vtkAreaSimplificationMetric* New();

  double ComputeMetric(vtkDataSet *mesh, vtkDataArray *scalarField,
    vtkIdType startCriticalPoint, vtkAbstractArray *vertexList,
    vtkIdType endCriticalPoint);

};

#endif
