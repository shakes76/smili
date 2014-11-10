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
#include "vtkDiffusionTensorGlyphFilter.h"

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "vtkObjectFactory.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkAppendPolyData.h"

vtkStandardNewMacro(vtkDiffusionTensorGlyphFilter);

vtkDiffusionTensorGlyphFilter::vtkDiffusionTensorGlyphFilter() : vtkProgrammableGlyphFilter()
{
  m_Resolution = 16;
}

vtkDiffusionTensorGlyphFilter::~vtkDiffusionTensorGlyphFilter()
{

}

double vtkDiffusionTensorGlyphFilter::computeAmplitude(std::vector<double> SH, double x, double y, double z, int lmax)
{
    double az = atan2(y, x);
    double el = acos(z);
    double val = 0.0;
    for (int l = 0; l <= lmax; l+=2)
    {
        // val += SH[vtkDiffusionTensorGlyphFilter::index(l,0)] * boost::math::legendre_p<double>(l, 0, z);
        val += SH[vtkDiffusionTensorGlyphFilter::index(l,0)] * boost::math::spherical_harmonic_r<double>(l, 0, el,az);
    }

    for (int m = 1; m <= lmax; m++)
    {
//        float caz = cos(m*az);
//        float saz = sin(m*az);
        for (int l = 2*((m+1)/2); l <= lmax; l+=2)
        {
        // float buf = boost::math::legendre_p<double>(l, 0, z);
        // val += SH[index(l,m)]*buf*caz;
        // val += SH[index(l,-m)]*buf*saz;
            val += SH[vtkDiffusionTensorGlyphFilter::index(l,m)]*boost::math::spherical_harmonic_r(l,m,el,az);
            val += SH[vtkDiffusionTensorGlyphFilter::index(l,-m)]*boost::math::spherical_harmonic_i(l,m,el,az);
        }
    }

    return val;
}

//Overide the vtkProgrammableGlyphFilter due to colouring not being unsigned char etc.
/*int vtkDiffusionTensorGlyphFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *source = vtkPolyData::SafeDownCast(
    sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPointData *inputPD = input->GetPointData();
  vtkCellData *inputCD = input->GetCellData();
  vtkPointData *outputPD = output->GetPointData();
  vtkCellData *outputCD = output->GetCellData();
  vtkPoints *newPts, *sourcePts;
  vtkUnsignedCharArray *ptScalars=NULL, *cellScalars=NULL; //CHANGED
  vtkDataArray *inPtScalars = NULL, *inCellScalars = NULL;
  vtkIdType numPts = input->GetNumberOfPoints();
  vtkPointData *sourcePD;
  vtkCellData *sourceCD;
  vtkIdType numSourcePts, numSourceCells, ptOffset=0, cellId, ptId, id, idx;
  int i, npts;
  vtkIdList *pts;
  vtkIdList *cellPts;
  vtkCell *cell;

  // Initialize
  vtkDebugMacro(<<"Generating programmable glyphs!");
  std::cout << "Generating Diffusion programmable glyphs!" << std::endl;

  if ( numPts < 1 )
    {
    vtkErrorMacro(<<"No input points to glyph");
    }

  pts=vtkIdList::New();
  pts->Allocate(VTK_CELL_SIZE);
  sourcePD = source->GetPointData();
  sourceCD = source->GetCellData();
  numSourcePts = source->GetNumberOfPoints();
  numSourceCells = source->GetNumberOfCells();

  outputPD->CopyScalarsOff(); //'cause we control the coloring process
  outputCD->CopyScalarsOff();

  output->Allocate(numSourceCells*numPts,numSourceCells*numPts);
  outputPD->CopyAllocate(sourcePD, numSourcePts*numPts, numSourcePts*numPts); //causes float array allocation
  outputCD->CopyAllocate(sourceCD, numSourceCells*numPts, numSourceCells*numPts);
  newPts = vtkPoints::New();
  newPts->Allocate(numSourcePts*numPts);

  // figure out how to color the data and setup
  if ( this->ColorMode == VTK_COLOR_BY_INPUT )
    {
    if ( (inPtScalars = inputPD->GetScalars()) )
      {
      ptScalars = vtkUnsignedCharArray::New();
      ptScalars->SetNumberOfComponents(3); //CHANGED
      ptScalars->SetName("Colours");
      ptScalars->SetNumberOfTuples(numSourcePts*numPts);
      }
    if ( (inCellScalars = inputCD->GetScalars()) )
      {
      cellScalars = vtkUnsignedCharArray::New();
      cellScalars->SetNumberOfComponents(3); //CHANGED
      cellScalars->SetNumberOfTuples(numSourcePts*numPts);
      }
    }

  else
    {
    if ( sourcePD->GetScalars() )
      {
      ptScalars = vtkUnsignedCharArray::New();
      ptScalars->SetNumberOfComponents(3); //CHANGED
      ptScalars->SetName("Colours");
      ptScalars->SetNumberOfTuples(numSourcePts*numPts);
      }
    if ( sourceCD->GetScalars() )
      {
      cellScalars = vtkUnsignedCharArray::New();
      cellScalars->SetNumberOfComponents(3); //CHANGED
      cellScalars->SetNumberOfTuples(numSourcePts*numPts);
      }
    }

  // Loop over all points, invoking glyph method and Update(),
  // then append output of source to output of this filter.
  //
//  this->Updating = 1; // to prevent infinite recursion
  this->PointData = input->GetPointData();
  for (this->PointId=0; this->PointId < numPts; this->PointId++)
    {
    if ( ! (this->PointId % 10000) )
      {
      this->UpdateProgress ((double)this->PointId/numPts);
      if (this->GetAbortExecute())
        {
        break;
        }
      }

    input->GetPoint(this->PointId, this->Point);

    if ( this->GlyphMethod )
      {
      (*this->GlyphMethod)(this->GlyphMethodArg);

      // The GlyphMethod may have set the source connection to NULL
      if (this->GetNumberOfInputConnections(1) == 0)
        {
        source = NULL;
        }
      else
        {
        // Update the source connection in case the GlyphMethod changed
        // its parameters.
        this->GetInputAlgorithm(1, 0)->Update();
        // The GlyphMethod may also have changed the source.
        sourceInfo = inputVector[1]->GetInformationObject(0);
        source = vtkPolyData::SafeDownCast(
          sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
        }
      }
    if (source)
      {
      sourcePts = source->GetPoints();
      numSourcePts = source->GetNumberOfPoints();
      numSourceCells = source->GetNumberOfCells();
      sourcePD = source->GetPointData();
      sourceCD = source->GetCellData();

      if ( this->ColorMode == VTK_COLOR_BY_SOURCE )
        {
        inPtScalars = sourcePD->GetScalars();
        inCellScalars = sourceCD->GetScalars();
        }

      // Copy all data from source to output.
      for (ptId=0; ptId < numSourcePts; ptId++)
        {
        id = newPts->InsertNextPoint(sourcePts->GetPoint(ptId));
        outputPD->CopyData(sourcePD, ptId, id);
        }

      for (cellId=0; cellId < numSourceCells; cellId++)
        {
        cell = source->GetCell(cellId);
        cellPts = cell->GetPointIds();
        npts = cellPts->GetNumberOfIds();
        for (pts->Reset(), i=0; i < npts; i++)
          {
          pts->InsertId(i,cellPts->GetId(i) + ptOffset);
          }
        id = output->InsertNextCell(cell->GetCellType(),pts);
        outputCD->CopyData(sourceCD, cellId, id);
        }

      // If we're coloring the output with scalars, do that now
      if ( ptScalars )
        {
        for (ptId=0; ptId < numSourcePts; ptId++)
          {
          idx = (this->ColorMode == VTK_COLOR_BY_INPUT ? this->PointId : ptId);
          unsigned char colourOfPoint[3] = {0, 0, 0};
          colourOfPoint[0] = static_cast<unsigned char>( inPtScalars->GetComponent(idx, 0) );
          colourOfPoint[1] = static_cast<unsigned char>( inPtScalars->GetComponent(idx, 1) );
          colourOfPoint[2] = static_cast<unsigned char>( inPtScalars->GetComponent(idx, 2) );
          ptScalars->InsertNextTupleValue(colourOfPoint); //CHANGED
          }
        }
      else if ( cellScalars )
        {
        for (cellId=0; cellId < numSourceCells; cellId++)
          {
          idx = (this->ColorMode == VTK_COLOR_BY_INPUT ? this->PointId : cellId);
          unsigned char colourOfCell[3] = {0, 0, 0};
          colourOfCell[0] = static_cast<unsigned char>( inCellScalars->GetComponent(idx, 0) );
          colourOfCell[1] = static_cast<unsigned char>( inCellScalars->GetComponent(idx, 1) );
          colourOfCell[2] = static_cast<unsigned char>( inCellScalars->GetComponent(idx, 2) );
          cellScalars->InsertNextTupleValue(colourOfCell); //CHANGED
          }
        }

      ptOffset += numSourcePts;

      }//if a source is available
    } //for all input points

//  this->Updating = 0;

  pts->Delete();

  output->SetPoints(newPts);
  newPts->Delete();

  if ( ptScalars )
    {
    idx = outputPD->AddArray(ptScalars);
    outputPD->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    ptScalars->Delete();
    }

  if ( cellScalars )
    {
    idx = outputCD->AddArray(cellScalars);
    outputCD->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    cellScalars->Delete();
    }

  output->Squeeze();

  return 1;
}*/
int vtkDiffusionTensorGlyphFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *source = vtkPolyData::SafeDownCast(
    sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPointData *inputPD = input->GetPointData();
  vtkCellData *inputCD = input->GetCellData();
  vtkIdType numPts = input->GetNumberOfPoints();
  vtkPointData *sourcePD;
  vtkCellData *sourceCD;
  vtkIdType numSourcePts, numSourceCells, ptOffset=0, cellId, ptId, id, idx;

  vtkSmartPointer<vtkAppendPolyData> appendedOutput = vtkSmartPointer<vtkAppendPolyData>::New();
  #if VTK_MAJOR_VERSION <= 5
    appendedOutput->AddInput(output);
  #else
    appendedOutput->AddInputData(output);
  #endif

  // Initialize
  vtkDebugMacro(<<"Generating programmable glyphs!");

  if ( numPts < 1 )
    {
    vtkErrorMacro(<<"No input points to glyph");
    }

  sourcePD = source->GetPointData();
  sourceCD = source->GetCellData();
  numSourcePts = source->GetNumberOfPoints();
  numSourceCells = source->GetNumberOfCells();

  // Loop over all points, invoking glyph method and Update(),
  // then append output of source to output of this filter.
  //
//  this->Updating = 1; // to prevent infinite recursion
  this->PointData = input->GetPointData();
  for (this->PointId=0; this->PointId < numPts; this->PointId++)
    {
    if ( ! (this->PointId % 10000) )
      {
      this->UpdateProgress (static_cast<long double>(this->PointId)/numPts);
      if (this->GetAbortExecute())
        {
        break;
        }
      }

    input->GetPoint(this->PointId, this->Point);

    if ( this->GlyphMethod )
      {
      (*this->GlyphMethod)(this->GlyphMethodArg);

    #if VTK_MAJOR_VERSION > 5
      // The GlyphMethod may have set the source connection to NULL
      if (this->GetNumberOfInputConnections(1) == 0)
        {
        source = NULL;
        }
      else
        {
        // Update the source connection in case the GlyphMethod changed
        // its parameters.
        this->GetInputAlgorithm(1, 0)->Update();
        // The GlyphMethod may also have changed the source.
        sourceInfo = inputVector[1]->GetInformationObject(0);
        source = vtkPolyData::SafeDownCast(
          sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
        }
    #endif
      }
    if (source)
      {
      numSourcePts = source->GetNumberOfPoints();
      numSourceCells = source->GetNumberOfCells();
      sourcePD = source->GetPointData();
      sourceCD = source->GetCellData();

      #if VTK_MAJOR_VERSION <= 5
        appendedOutput->AddInput(source);
      #else
        appendedOutput->AddInputData(source);
      #endif

      ptOffset += numSourcePts;
      }//if a source is available
    } //for all input points

//  this->Updating = 0;
  appendedOutput->Update();
  this->UpdateProgress(100);
  output->DeepCopy(appendedOutput->GetOutput());
  output->Squeeze();

  return 1;
}

void vtkDiffusionTensorGlyphFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Color Mode: " << this->GetColorModeAsString() << endl;
  os << indent << "Point Id: " << this->PointId << "\n";
  os << indent << "Point: " << this->Point[0]
                    << ", " << this->Point[1]
                    << ", " << this->Point[2] << "\n";
  if (this->PointData)
    {
    os << indent << "PointData: " << this->PointData << "\n";
    }
  else
    {
    os << indent << "PointData: (not defined)\n";
    }

  if ( this->GlyphMethod )
    {
    os << indent << "Glyph Method defined\n";
    }
  else
    {
    os << indent << "No Glyph Method\n";
    }
}
