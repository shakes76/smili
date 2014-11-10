/*=========================================================================
  Program: milxSMILI
  Module: milxDeformableModel
  Author: Shekhar Chandra
  Language: C++
  Created: 18 Apr 2011 12:03:00 EST

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
#include "milxDeformableModel.h"
//SMILI
#include "milxMath.h"

namespace milx
{

DeformableModel::DeformableModel() : Model()
{
  InternalInPlaceOperation = false;
}

DeformableModel::DeformableModel(vtkSmartPointer<vtkPolyData> model)
{
  InternalInPlaceOperation = false;
  Model::SetInput(model);
}

} //end milx namespace
